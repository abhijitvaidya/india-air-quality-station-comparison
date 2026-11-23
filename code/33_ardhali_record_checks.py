#!/usr/bin/env python3
"""Checks on the Ardhali Bazar record, from hourly data 2016 to September 2025 (Section 4.5, Appendix G, Tables G1 to G3):
daily cycle by year; correlation of daily anomalies with three Bihar stations; response to rain; missing hours; ratios to
Patna for each pollutant in the twelve months before and after July 2021; yearly medians of the gas channels.

Run:  python code/33_ardhali_record_checks.py   (needs data/raw/aq.duckdb and data/raw/weather; see data/README.md)
Outputs: results/ardhali_record_checks.md and results/ardhali_*.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os
import numpy as np, pandas as pd, duckdb

HERE = str(RESULTS)
ARD = 'site_273'
NEW = {'site_5469': 'Bhelupur', 'site_5468': 'IESD BHU', 'site_5475': 'Maldahiya'}
BIH = {'site_157': 'Patna IGSC', 'site_159': 'Gaya', 'site_158': 'Muzaffarpur'}
UPC = {'site_277': 'Lucknow Lalbagh', 'site_297': 'Lucknow Talkatora', 'site_276': 'Kanpur Nehru Nagar'}
NAMES = {ARD: 'Ardhali Bazar', **NEW, **BIH, **UPC}
CITY = {ARD: 'Varanasi', **{k: 'Varanasi' for k in NEW}, 'site_157': 'Patna', 'site_159': 'Patna', 'site_158': 'Patna',
        'site_277': 'Lucknow', 'site_297': 'Lucknow', 'site_276': 'Kanpur'}
EXTRA = {'site_5262': 'Patna Rajbansi Nagar', 'site_5261': 'Patna Muradpur'}   # second Patna stations for intra-city pairs
DIWALI = {2016: '2016-10-30', 2017: '2017-10-19', 2018: '2018-11-07', 2019: '2019-10-27', 2020: '2020-11-14',
          2021: '2021-11-04', 2022: '2022-10-24', 2023: '2023-11-12', 2024: '2024-10-31'}

c = duckdb.connect(str(HOURLY_DB), read_only=True)
ids = list(NAMES) + list(EXTRA)
h = c.sql(f"""select station_id sid, p, date_trunc('hour', ts) ts, avg(v) v from pm
  where source = 'cpcb_caaqm' and station_id in {tuple(ids)} and ts >= '2016-01-01' and v > 0 and v <= 1000
  group by all""").df()
raw = c.sql(f"""select station_id sid, v from pm where source = 'cpcb_caaqm' and p = 'PM2.5'
  and station_id in {tuple(ids)} and ts >= '2016-01-01' and v > 0 and v < 985""").df()
h['y'] = h.ts.dt.year; h['d'] = h.ts.dt.normalize(); h['hr'] = h.ts.dt.hour
pm25 = h[h.p == 'PM2.5'].copy()
L = ['# Checks on the Ardhali Bazar record, before and after 2021', '',
     'Hourly CPCB PM2.5 (XKDR), 2016 to Sep 2025. Ardhali Bazar against three new Varanasi stations (from Jun-Jul 2021), '
     'three long-record Bihar stations and three long-record UP stations. These are checks on the measurements only.', '']


def table(df, title, fmt='{:.2f}'):
    global L
    L += [f'## {title}', '', '| station | ' + ' | '.join(str(c) for c in df.columns) + ' |', '|---|' + '---|' * len(df.columns)]
    for idx, r in df.iterrows():
        L.append(f'| {NAMES.get(idx, EXTRA.get(idx, idx))} | ' + ' | '.join('' if pd.isna(v) else fmt.format(v) for v in r.values) + ' |')
    L.append('')


order = [ARD] + list(NEW) + list(BIH) + list(UPC)

# 1 diurnal cycle ------------------------------------------------------------------------------
full = pm25.groupby(['sid', 'd']).v.transform('count') >= 20
p1 = pm25[full].copy()
p1['rel'] = p1.v / p1.groupby(['sid', 'd']).v.transform('mean')
prof = p1.groupby(['sid', 'y', 'hr']).rel.mean().unstack('hr')
amp = (prof.max(axis=1) / prof.min(axis=1)).unstack('y')
base = p1[p1.y.between(2016, 2018)].groupby(['sid', 'hr']).rel.mean().unstack('hr')
shape = pd.DataFrame({(s, y): [np.corrcoef(prof.loc[(s, y)], base.loc[s])[0, 1]] for s, y in prof.index if s in base.index}).T
shape = shape[0].unstack(1) if len(shape) else shape
peak = prof.idxmax(axis=1).unstack('y'); trough = prof.idxmin(axis=1).unstack('y')
prof.to_csv(f'{HERE}/ardhali_diurnal_profiles.csv')
table(amp.reindex(order), '1a. Day-night amplitude (highest hour / lowest hour of the mean daily profile)')
table(shape.reindex(order), "1b. Correlation of each year's daily profile with the station's own 2016-18 profile")
table(peak.reindex(order).astype('Int64'), '1c. Hour of the daily peak', '{}')
table(trough.reindex(order).astype('Int64'), '1d. Hour of the daily low', '{}')

# 2 day-to-day coupling -------------------------------------------------------------------------
daily = pm25.groupby(['sid', 'd']).v.agg(['mean', 'count']).reset_index()
daily = daily[daily['count'] >= 16]
daily['lv'] = np.log(daily['mean'])
daily = daily.sort_values(['sid', 'd'])
daily['anom'] = daily.groupby('sid').lv.transform(lambda s: s - s.rolling(31, center=True, min_periods=15).median())
A = daily.pivot(index='d', columns='sid', values='anom')
A.to_csv(f'{HERE}/ardhali_daily_anomalies.csv')
rows = {}
for s in order:
    if s not in A: continue
    refcols = [b for b in BIH if b != s and b in A]
    ref = A[refcols].median(axis=1)
    ok = A[s].notna() & ref.notna()
    rows[s] = {y: A[s][ok & (A.index.year == y)].corr(ref[ok & (A.index.year == y)]) for y in range(2016, 2026)
               if (ok & (A.index.year == y)).sum() >= 60}
table(pd.DataFrame(rows).T.reindex(order), '2. Day-to-day coupling: correlation of daily anomalies with the Bihar reference (Patna, Gaya, Muzaffarpur; median, excluding the station itself)')

# 3 rain response -------------------------------------------------------------------------------
w = pd.read_csv(WEATHER_CSV, parse_dates=['d'])
rr = {}
for s in order:
    if s not in A: continue
    ww = w[w.city == CITY[s]].set_index('d').rain_mm
    x = pd.DataFrame({'a': A[s], 'rain': ww}).dropna()
    out = {}
    for lab, yrs in [('2016-19', range(2016, 2020)), ('2022-24', range(2022, 2025))]:
        xx = x[x.index.year.isin(list(yrs))]
        wet, dry = xx[xx.rain >= 10].a, xx[xx.rain < 1].a
        out[f'{lab} wet-dry'] = (np.exp(wet.mean() - dry.mean()) - 1) * 100 if len(wet) >= 10 else np.nan
        out[f'{lab} wet days'] = len(wet)
    rr[s] = out
table(pd.DataFrame(rr).T.reindex(order), "3. Rain response: % difference in the daily anomaly on days with >= 10 mm rain in the station's city vs dry days", '{:.0f}')

# 4 Diwali --------------------------------------------------------------------------------------
dv = {}
for s in order:
    x = pm25[pm25.sid == s].set_index('ts').v
    out = {}
    for y, ds in DIWALI.items():
        d0 = pd.Timestamp(ds)
        def night(dd):
            seg = x[(x.index >= dd + pd.Timedelta(hours=18)) & (x.index < dd + pd.Timedelta(hours=30))]
            return seg.mean() if len(seg) >= 6 else np.nan
        base_n = np.nanmean([night(d0 - pd.Timedelta(days=k)) for k in range(3, 8)])
        out[y] = night(d0) / base_n if base_n and not np.isnan(base_n) else np.nan
    dv[s] = out
table(pd.DataFrame(dv).T.reindex(order), '4. Diwali night (18:00-05:59) over the mean of nights D-7..D-3')

# 5 texture -------------------------------------------------------------------------------------
tx = {}
for s in order:
    x = pm25[pm25.sid == s].sort_values('ts')
    out = {}
    for y in range(2016, 2026):
        xy = x[x.y == y]
        if len(xy) < 2000: continue
        lv = np.log(xy.set_index('ts').v)
        an = lv - lv.groupby(lv.index.normalize()).transform('mean')
        cont = an.index.to_series().diff() == pd.Timedelta(hours=1)
        out[(y, 'sd_daily_anom')] = A[s][A.index.year == y].std() if s in A else np.nan
        out[(y, 'lag1_hourly')] = an[cont].corr(an.shift(1)[cont])
        out[(y, 'repeat_share')] = ((xy.v.diff() == 0) & cont.values).mean()
    tx[s] = out
T = pd.DataFrame(tx).T
for k in ['sd_daily_anom', 'lag1_hourly', 'repeat_share']:
    sub = T.xs(k, axis=1, level=1) if isinstance(T.columns, pd.MultiIndex) else None
    if sub is not None: table(sub.reindex(order), f'5. Texture: {k}', '{:.3f}')
# digits and whole numbers from raw values
raw['y'] = 0
dg = {}
for s in order:
    v = raw[raw.sid == s].v.values
    if len(v) < 1000: continue
    whole = np.mean(np.isclose(v, np.round(v)))
    ld = (np.floor(np.round(v * 100)) % 10).astype(int)          # second decimal digit
    li = (np.floor(v) % 10).astype(int)                           # last integer digit
    def chi(d):
        obs = np.bincount(d, minlength=10); e = obs.sum() / 10
        return ((obs - e) ** 2 / e).sum()
    dg[s] = {'n values': len(v), 'whole-number share': whole, 'chi2 last integer digit (9 df; 5% crit 16.9)': chi(li),
             'chi2 2nd decimal digit': chi(ld)}
table(pd.DataFrame(dg).T.reindex(order), '5b. Raw-value digits, 2016-2025 (PM2.5)', '{:.3f}')
# PM2.5 / PM10
pp = h.pivot_table(index=['sid', 'ts'], columns='p', values='v').dropna().reset_index()
pp['y'] = pp.ts.dt.year; pp['r'] = pp['PM2.5'] / pp['PM10']
pr = pp.groupby(['sid', 'y']).agg(ratio_mean=('r', 'mean'), ratio_sd=('r', 'std'), n=('r', 'size'))
pr['logcorr'] = pp.groupby(['sid', 'y']).apply(lambda g: np.corrcoef(np.log(g['PM2.5']), np.log(g['PM10']))[0, 1])
for k in ['ratio_mean', 'ratio_sd', 'logcorr']:
    table(pr[k].unstack('y').reindex(order), f'5c. PM2.5/PM10 hourly: {k}', '{:.2f}')

# 6 co-movement after July 2021 -----------------------------------------------------------------
hp = pm25.pivot_table(index='ts', columns='sid', values='v')
hp = hp[hp.index >= '2021-08-01']
la = np.log(hp)
la = la - la.groupby(la.index.normalize()).transform('mean')     # within-day hourly anomaly
pairs = [(ARD, k) for k in NEW] + [('site_5469', 'site_5468'), ('site_5469', 'site_5475'), ('site_5468', 'site_5475'),
         ('site_277', 'site_297'), ('site_157', 'site_5262'), ('site_157', 'site_5261')]
cm = []
for a, b in pairs:
    if a in la and b in la:
        ok = la[a].notna() & la[b].notna()
        for y in range(2021, 2026):
            m = ok & (la.index.year == y)
            if m.sum() > 500:
                cm.append({'pair': f'{NAMES.get(a, EXTRA.get(a))} / {NAMES.get(b, EXTRA.get(b))}', 'year': y,
                           'hourly r': la[a][m].corr(la[b][m]), 'daily r': A[a].corr(A[b]) if a in A and b in A else np.nan,
                           'level ratio (median)': (hp[a][m] / hp[b][m]).median()})
cmd = pd.DataFrame(cm)
cmd.to_csv(f'{HERE}/ardhali_comovement.csv', index=False)
L += ['## 6. Co-movement from Aug 2021: within-day hourly anomaly correlation and level ratio', '',
      '| pair | year | hourly r | level ratio a/b (median) |', '|---|---|---|---|']
L += [f'| {r.pair} | {r.year} | {r["hourly r"]:.2f} | {r["level ratio (median)"]:.2f} |' for _, r in cmd.iterrows()] + ['']

# 7 missing hours -------------------------------------------------------------------------------
ms = {}
for s in order:
    x = pm25[pm25.sid == s]
    out = {}
    for y in range(2016, 2026):
        hrs = 8760 if y < 2025 else 5856
        out[y] = 1 - x[x.y == y].ts.nunique() / hrs
    ms[s] = out
table(pd.DataFrame(ms).T.reindex(order), '7a. Share of hours missing (PM2.5)', '{:.2f}')
# are Ardhali's missing hours the high ones? compare new-station level when Ardhali missing vs present (from Aug 2021)
nm = hp[[k for k in NEW if k in hp]].median(axis=1)
miss = hp[ARD].isna()
L += [f'7b. From Aug 2021, median PM2.5 at the three new Varanasi stations: when Ardhali is missing '
      f'{nm[miss].median():.0f}, when present {nm[~miss].median():.0f} (hours: {int(miss.sum())} missing, {int((~miss).sum())} present).', '']

# 8 gases at Ardhali ---------------------------------------------------------------------------
g = c.sql(f"""select station_id sid, parameter_name p, date_trunc('hour', collected_at) ts, avg(value) v from read_parquet('{RAW}/xkdr/measurements/**/*.parquet', hive_partitioning=true)
  where station_id in ('{ARD}', 'site_157') and parameter_name in ('NO2', 'NOx', 'NO', 'CO', 'SO2', 'Ozone', 'PM2.5', 'PM10')
    and collected_at >= '2016-01-01' and value > 0 group by all""").df()
g['y'] = g.ts.dt.year
med = g.groupby(['sid', 'p', 'y']).v.median().unstack('y')
med.to_csv(f'{HERE}/ardhali_gas_annual_medians.csv')
L += ['## 8a. Annual medians, Ardhali Bazar (hourly; µg/m³, CO mg/m³, NOx ppb)', '',
      '| pollutant | ' + ' | '.join(str(y) for y in med.columns) + ' |', '|---|' + '---|' * len(med.columns)]
for p in ['PM2.5', 'PM10', 'NO2', 'NO', 'NOx', 'CO', 'SO2', 'Ozone']:
    if (ARD, p) in med.index:
        L.append(f'| {p} | ' + ' | '.join('' if pd.isna(v) else f'{v:.1f}' for v in med.loc[(ARD, p)].values) + ' |')
L.append('')
ng = g[g.sid == ARD].pivot_table(index='ts', columns='p', values='v')
if {'NO2', 'NOx'} <= set(ng.columns):
    x = ng[['NO2', 'NOx']].dropna()
    bad = (x.NO2 / 1.88 > x.NOx * 1.02)
    L += ['8b. Share of hours with NO2 (converted to ppb) above NOx (ppb), Ardhali: ' +
          ', '.join(f'{y}: {bad[x.index.year == y].mean():.0%}' for y in range(2016, 2026) if (x.index.year == y).sum() > 500), '']
# monthly ratios to Patna IGSC around July 2021
mo = g.assign(m=g.ts.dt.to_period('M')).groupby(['sid', 'p', 'm']).v.median().unstack('sid').dropna()
mo['ratio'] = mo[ARD] / mo['site_157']
rt = mo['ratio'].unstack('p')
win = rt[(rt.index >= pd.Period('2020-07', 'M')) & (rt.index <= pd.Period('2022-06', 'M'))]
win.to_csv(f'{HERE}/ardhali_ratio_to_patna_monthly.csv')
pre = rt[(rt.index >= pd.Period('2020-07', 'M')) & (rt.index <= pd.Period('2021-06', 'M'))].median()
post = rt[(rt.index >= pd.Period('2021-08', 'M')) & (rt.index <= pd.Period('2022-07', 'M'))].median()
L += ['## 8c. Ardhali / Patna IGSC monthly median ratio, 12 months before (Jul 2020-Jun 2021) vs after (Aug 2021-Jul 2022) the July 2021 step', '',
      '| pollutant | before | after | change |', '|---|---|---|---|']
for p in pre.index:
    if p in post.index and pd.notna(pre[p]) and pd.notna(post[p]):
        L.append(f'| {p} | {pre[p]:.2f} | {post[p]:.2f} | {100 * (post[p] / pre[p] - 1):+.0f}% |')
L.append('')
open(f'{HERE}/ardhali_record_checks.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
