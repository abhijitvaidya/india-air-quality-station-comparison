#!/usr/bin/env python3
"""Monthly series, linked year-on-year check and steady-ratio test (Figure 1, Table 3, Table E3).

1. Monthly series: the calculation of 11_national_series.py run one calendar month at a time, so each month is compared
   with the same month of 2016-2019; a centred 12-month mean is added (missing within six months of either end).
2. Linked check: for each pair of adjacent years, stations with at least six common calendar months; their mean change
   in log ratio is chained from year to year and compared with the main estimates.
3. Steady-ratio test: each station's mean change in log ratio between two years over the calendar months present in
   both (at least six); 2019 to 2022 for stations reporting by 2019, and 2022 to 2024 for those stations and for
   stations that first report in 2020 or later. Intervals resample cities.

Run after 11_national_series.py:  python code/12_national_monthly.py
Outputs: results/national_monthly.csv/.md, results/national_matched.csv, results/national_steady.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys
import numpy as np, pandas as pd
import netCDF4 as nc

H = str(RESULTS)
from common import paired

YEARS = list(range(2016, 2025)); BASE = [2016, 2017, 2018, 2019]; REF = 2019; NB = 2000
rng = np.random.default_rng(20261122)
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']
RN = {'washu': 'WashU', 'cams': 'CAMS', 'ceres': 'CERES'}

m = paired('v6')
st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id')
d = nc.Dataset(str(CAMS_NC))
tv = nc.num2date(d['valid_time'][:], d['valid_time'].units, 'proleptic_gregorian')
YM = [f'{t.year}{t.month:02d}' for t in tv]; LAT = np.array(d['latitude'][:]); LON = np.array(d['longitude'][:]); F = np.array(d['pm2p5'][:]) * 1e9
la_ = LAT if LAT[0] < LAT[-1] else LAT[::-1]; f_ = F if LAT[0] < LAT[-1] else F[:, ::-1, :]
def bilinear(lat, lon):
    i = np.clip(np.searchsorted(la_, lat) - 1, 0, len(la_) - 2); j = np.clip(np.searchsorted(LON, lon) - 1, 0, len(LON) - 2)
    wy = (lat - la_[i]) / (la_[i + 1] - la_[i]); wx = (lon - LON[j]) / (LON[j + 1] - LON[j])
    return (1 - wy) * (1 - wx) * f_[:, i, j] + (1 - wy) * wx * f_[:, i, j + 1] + wy * (1 - wx) * f_[:, i + 1, j] + wy * wx * f_[:, i + 1, j + 1]
cams, aod = {}, {}
for sid in m.station_id.unique():
    la, lo = st.lat[sid], st.lon[sid]
    for ym, v in zip(YM, bilinear(la, lo)): cams[(sid, ym)] = v
    f = aod_file(la, lo)
    if os.path.exists(f):
        for k, v in json.load(open(f))['properties']['parameter']['AOD_55'].items():
            if k[4:] != '13' and v is not None and v > 0: aod[(sid, k)] = v
m['l_st'] = np.log(m.cpcb); m['l_washu'] = np.log(m.sat_pm25)
m['l_cams'] = np.log(pd.Series([cams.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float))
m['l_ceres'] = np.log(pd.Series([aod.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float))
m['citykey'] = m.state + '/' + m.city.fillna(m.station_id)
up = m.state == 'Uttar Pradesh'
first = m.groupby('station_id').y.min(); m['new'] = m.station_id.map(first) >= 2020
GROUPS = {'Uttar Pradesh': up, 'Rest of India': ~up, 'Plain outside Uttar Pradesh': m.state.isin(IGP)}
bi = [YEARS.index(y) for y in BASE]; keep = [i for i, y in enumerate(YEARS) if y != REF]


def solve(As, Bs, nref):
    dead = np.einsum('...ii->...i', As) < 1e-9
    beta = np.where(dead, np.nan, np.einsum('...ij,...j->...i', np.linalg.pinv(As, hermitian=True), Bs))
    full = np.insert(beta, YEARS.index(REF), 0.0, axis=-1)
    full = full - np.nanmean(full[..., bi], axis=-1, keepdims=True)
    return np.where(np.asarray(nref)[..., None] == 0, np.nan, full)


def monthly(g, ycol):
    """Returns est[108], draws[NB,108], stations by month[108]. Same city draws for all calendar months."""
    g = g.dropna(subset=[ycol]); keys, inv = np.unique(g.citykey.values, return_inverse=True)
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    est = np.full(108, np.nan); dr = np.full((NB, 108), np.nan); cnt = np.zeros(108, int)
    for mo in range(1, 13):
        sel = g.mo.values == mo; s = g[sel]; iv = inv[sel]
        X = np.zeros((len(s), len(YEARS))); X[np.arange(len(s)), s.y.values - YEARS[0]] = 1.0
        yv = s[ycol].values - s.groupby('station_id')[ycol].transform('mean').values
        Xd = (X - pd.DataFrame(X).groupby(s.station_id.values).transform('mean').values)[:, keep]
        A = np.zeros((len(keys), len(keep), len(keep))); B = np.zeros((len(keys), len(keep))); nref = np.zeros(len(keys))
        for k in np.unique(iv):
            idx = iv == k; A[k] = Xd[idx].T @ Xd[idx]; B[k] = Xd[idx].T @ yv[idx]; nref[k] = (s.y.values[idx] == REF).sum()
        ix = [(y - YEARS[0]) * 12 + mo - 1 for y in YEARS]
        est[ix] = solve(A.sum(0), B.sum(0), nref.sum()); dr[:, ix] = solve(np.tensordot(W, A, 1), W @ B, W @ nref)
        cnt[ix] = s.groupby('y').station_id.nunique().reindex(YEARS).fillna(0).astype(int).values
    return est, dr, cnt


def smooth(a):
    """Centred 12-month mean (2 x 12); NaN within six months of either end or where a month is missing."""
    k = np.r_[0.5, np.ones(11), 0.5] / 12
    return np.apply_along_axis(lambda v: pd.Series(v).rolling(13, center=True).apply(lambda w: np.dot(w, k), raw=True).values, -1, a)


pc = lambda v: 100 * np.expm1(v)
rows = []
for gname, mask in GROUPS.items():
    g = m[mask]
    for ref in ['washu', 'cams', 'ceres']:
        gg = g.dropna(subset=[f'l_{ref}']).copy(); gg['l_ratio'] = gg.l_st - gg[f'l_{ref}']
        for what, col in ([('station', 'l_st'), ('reference', f'l_{ref}'), ('ratio', 'l_ratio')] if ref == 'washu' else [('ratio', 'l_ratio')]):
            est, dr, cnt = monthly(gg, col); es = smooth(est); ds = smooth(dr)
            for t in range(108):
                rows.append({'group': gname, 'reference': ref, 'series': what, 'ym': f'{2016 + t // 12}{t % 12 + 1:02d}', 'pct': pc(est[t]),
                             'lo': pc(np.nanpercentile(dr[:, t], 2.5)) if np.isfinite(dr[:, t]).any() else np.nan,
                             'hi': pc(np.nanpercentile(dr[:, t], 97.5)) if np.isfinite(dr[:, t]).any() else np.nan,
                             'pct_12m': pc(es[t]), 'lo_12m': pc(np.nanpercentile(ds[:, t], 2.5)) if np.isfinite(ds[:, t]).any() else np.nan,
                             'hi_12m': pc(np.nanpercentile(ds[:, t], 97.5)) if np.isfinite(ds[:, t]).any() else np.nan, 'stations': int(cnt[t])})
    print('monthly', gname, 'done', flush=True)
M = pd.DataFrame(rows); M.round(2).to_csv(f'{H}/national_monthly.csv', index=False)


# ---------- matched changes ----------
def matched(g, col, y0, y1, min_months=6):
    """Per-station mean change in col between y0 and y1 over common calendar months. Returns Series indexed by station."""
    a = g[g.y == y0].set_index(['station_id', 'mo'])[col]; b = g[g.y == y1].set_index(['station_id', 'mo'])[col]
    dd = (b - a).dropna().reset_index()
    s = dd.groupby('station_id')[col].agg(['mean', 'size'])
    return s[s['size'] >= min_months]['mean']


def group_mean(s, citymap, keys, W):
    """Mean over stations and its city-bootstrap draws (shared W over keys)."""
    if s.empty: return np.nan, np.full(len(W), np.nan), 0
    ck = citymap[s.index]; sm = pd.Series(s.values, index=ck.values).groupby(level=0).agg(['sum', 'size']).reindex(keys).fillna(0)
    num = W @ sm['sum'].values; den = W @ sm['size'].values
    return s.mean(), np.where(den > 0, num / np.where(den > 0, den, 1), np.nan), len(s)


citymap = m.groupby('station_id').citykey.first()
FE = pd.read_csv(f'{H}/national_series.csv'); FE = FE[FE.weights == 'station-months']
rows = []
for gname, mask in GROUPS.items():
    g0 = m[mask]; keys = np.unique(g0.citykey.values)
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    for ref in ['washu', 'cams', 'ceres']:
        g = g0.dropna(subset=[f'l_{ref}']).copy(); g['l_ratio'] = g.l_st - g[f'l_{ref}']
        steps, sdr, ns = [0.0], [np.zeros(NB)], [0]
        for y in YEARS[1:]:
            e, dr, n = group_mean(matched(g, 'l_ratio', y - 1, y), citymap, keys, W); steps.append(e); sdr.append(dr); ns.append(n)
        cum = np.cumsum(steps); cdr = np.cumsum(np.array(sdr), axis=0)          # [year, draw]
        cum = cum - cum[bi].mean(); cdr = cdr - cdr[bi].mean(0)
        fe = FE[(FE.group == gname) & (FE.reference == ref) & (FE.series == 'ratio')].set_index('year').pct
        for i, y in enumerate(YEARS):
            rows.append({'group': gname, 'reference': ref, 'year': y, 'stations_in_pair': ns[i], 'step_pct': pc(steps[i]), 'chained_pct': pc(cum[i]),
                         'lo': pc(np.nanpercentile(cdr[i], 2.5)), 'hi': pc(np.nanpercentile(cdr[i], 97.5)), 'main_series_pct': fe.get(y, np.nan)})
MT = pd.DataFrame(rows); MT.round(2).to_csv(f'{H}/national_matched.csv', index=False)

# ---------- steady-ratio test ----------
rows = []
SG = {'Uttar Pradesh, reporting by 2019': up & ~m.new, 'Uttar Pradesh, started 2020 or later': up & m.new,
      'Rest of India, reporting by 2019': ~up & ~m.new, 'Rest of India, started 2020 or later': ~up & m.new}
for gname, mask in SG.items():
    g0 = m[mask]; keys = np.unique(g0.citykey.values)
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    for ref in ['washu', 'cams', 'ceres']:
        g = g0.dropna(subset=[f'l_{ref}']).copy(); g['l_ratio'] = g.l_st - g[f'l_{ref}']
        for y0, y1 in [(2019, 2022), (2022, 2024)] + [(y - 1, y) for y in YEARS[1:]]:
            s = matched(g, 'l_ratio', y0, y1); e, dr, n = group_mean(s, citymap, keys, W)
            if n == 0: continue
            rows.append({'group': gname, 'reference': ref, 'from': y0, 'to': y1, 'stations': n, 'cities': citymap[s.index].nunique(), 'mean_pct': pc(e),
                         'lo': pc(np.nanpercentile(dr, 2.5)), 'hi': pc(np.nanpercentile(dr, 97.5)), 'median_pct': pc(s.median()), 'share_falling': (s < 0).mean()})
S = pd.DataFrame(rows); S.round(3).to_csv(f'{H}/national_steady.csv', index=False)
lv = m[m.y >= 2023].groupby('station_id').agg(lr=('lr', 'mean'), n=('lr', 'size')); lv = lv[lv.n >= 12]; lv['r'] = np.exp(lv.lr)
lv['grp'] = [next(k for k, v in SG.items() if sid in set(m[v].station_id)) for sid in lv.index]

# ---------- report ----------
fm = lambda v: 'n/a' if not np.isfinite(v) else f'{v:+.0f}%'.replace('-', '−')
L = ['# Monthly series, matched check and steady-ratio test', '',
     f'{len(m):,} station-months, {m.station_id.nunique()} stations. {int((first >= 2020).sum())} stations started in 2020 or later, {int(m[up & m.new].station_id.nunique())} of them in Uttar Pradesh.', '',
     '## 1. Monthly series: station against WashU, 12-month centred mean (selected months)', '', '| Month | Uttar Pradesh | Rest of India | Plain outside Uttar Pradesh |', '|---|---|---|---|']
def mcell(g, ym):
    x = M[(M.group == g) & (M.reference == 'washu') & (M.series == 'ratio') & (M.ym == ym)].iloc[0]
    return f'{fm(x.pct_12m)} ({fm(x.lo_12m)} to {fm(x.hi_12m)})'
for ym in ['201612', '201712', '201812', '201906', '201912', '202006', '202012', '202106', '202112', '202206', '202212', '202306', '202312', '202406']:
    L.append(f'| {ym[:4]}-{ym[4:]} | ' + ' | '.join(mcell(g, ym) for g in GROUPS) + ' |')
L += ['', 'Yearly mean of the monthly values against the main yearly series (station against WashU):', '', '| Year | UP, monthly method | UP, main series | Rest, monthly method | Rest, main series |', '|---|---|---|---|---|']
for y in YEARS:
    c = []
    for g in ['Uttar Pradesh', 'Rest of India']:
        x = M[(M.group == g) & (M.reference == 'washu') & (M.series == 'ratio') & (M.ym.str[:4] == str(y))]
        c += [fm(pc(np.log1p(x.pct / 100).mean())), fm(FE[(FE.group == g) & (FE.reference == 'washu') & (FE.series == 'ratio') & (FE.year == y)].pct.iloc[0])]
    L.append(f'| {y} | ' + ' | '.join(c) + ' |')
L += ['', '## 2. Matched year-on-year check: station against reference, chained, relative to the 2016–2019 mean', '',
      'Each step uses only stations with at least six common calendar months in both years. "Main" is the series of national_series.csv.', '']
for ref in ['washu', 'cams', 'ceres']:
    L += [f'### Against {RN[ref]}', '', '| Year | UP: stations in pair | UP chained | UP main | Rest: stations in pair | Rest chained | Rest main | Plain chained | Plain main |', '|---|---|---|---|---|---|---|---|---|']
    for y in YEARS:
        r = lambda g: MT[(MT.group == g) & (MT.reference == ref) & (MT.year == y)].iloc[0]
        u, o, p = r('Uttar Pradesh'), r('Rest of India'), r('Plain outside Uttar Pradesh')
        L.append(f'| {y} | {int(u.stations_in_pair)} | {fm(u.chained_pct)} ({fm(u.lo)} to {fm(u.hi)}) | {fm(u.main_series_pct)} | {int(o.stations_in_pair)} | {fm(o.chained_pct)} ({fm(o.lo)} to {fm(o.hi)}) | {fm(o.main_series_pct)} | {fm(p.chained_pct)} | {fm(p.main_series_pct)} |')
    L.append('')
L += ['## 3. Steady-ratio test: change in each station\'s ratio to the reference, same station and same calendar months', '',
      '| Group | Period | Reference | Stations (cities) | Mean change | Median | Share of stations falling |', '|---|---|---|---|---|---|---|']
for gname in SG:
    for y0, y1 in [(2019, 2022), (2022, 2024)]:
        for ref in ['washu', 'cams', 'ceres']:
            x = S[(S.group == gname) & (S.reference == ref) & (S['from'] == y0) & (S.to == y1)]
            if x.empty: continue
            x = x.iloc[0]; L.append(f'| {gname} | {y0} to {y1} | {RN[ref]} | {int(x.stations)} ({int(x.cities)}) | {fm(x.mean_pct)} ({fm(x.lo)} to {fm(x.hi)}) | {fm(x.median_pct)} | {x.share_falling:.0%} |')
L += ['', 'Year-on-year, against WashU (mean change; stations):', '', '| Year | ' + ' | '.join(SG) + ' |', '|---|' + '---|' * len(SG)]
for y in YEARS[1:]:
    c = []
    for gname in SG:
        x = S[(S.group == gname) & (S.reference == 'washu') & (S['from'] == y - 1) & (S.to == y)]
        c.append('' if x.empty else f'{fm(x.mean_pct.iloc[-1])} ({int(x.stations.iloc[-1])})')
    L.append(f'| {y} | ' + ' | '.join(c) + ' |')
L += ['', '## Level in 2023–2024: station PM2.5 as a multiple of WashU (stations with at least 12 valid months)', '', '| Group | Stations | 10th | 25th | Median | 75th | 90th | Share below 0.75 |', '|---|---|---|---|---|---|---|---|']
for gname in SG:
    q = lv[lv.grp == gname].r
    L.append(f'| {gname} | {len(q)} | ' + ' | '.join(f'{v:.2f}' for v in q.quantile([.1, .25, .5, .75, .9])) + f' | {(q < 0.75).mean():.0%} |')
open(f'{H}/national_monthly.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
