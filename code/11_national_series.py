#!/usr/bin/env python3
"""Year-by-year comparison of continuous-station PM2.5 with three references, 2016-2024 (Section 3.1, Table 2, Table E1).

Every valid station-month (at least 15 valid days) from 2016 to 2024 is used. A station is compared only with itself:
each station x calendar-month cell has its own mean, so a station that reads high or low for its surroundings does not
move the result, and a station that starts in 2021 informs the years from 2021 only. Year effects are estimated by least
squares on the cell-demeaned data, separately for each group of stations, relative to the mean of the 2016-2019
effects. The same estimator is applied to log(station PM2.5), log(reference at the station) and log(station / reference)
on identical station-months. Intervals: city bootstrap (2,000 draws), percentile 95%.

Groups: Uttar Pradesh; rest of India; rest of India without Delhi; Indo-Gangetic plain outside Uttar Pradesh (Punjab,
Haryana, Chandigarh, Delhi, Bihar, West Bengal). Variant: cities weighted equally.
References: satellite PM2.5 (WashU V6.GL.03), CAMS EAC4 PM2.5 (bilinear), CERES SYN1deg aerosol optical depth (NASA POWER).

Run: python code/11_national_series.py   (about 15 minutes)
Outputs: results/national_series.csv, results/national_series.md, results/national_levels.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys
import numpy as np, pandas as pd
import netCDF4 as nc

H = str(RESULTS)
from common import paired

YEARS = list(range(2016, 2025)); BASE = [2016, 2017, 2018, 2019]; NB = 2000
rng = np.random.default_rng(20261120)
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']

m = paired('v6')
st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id')
# CAMS PM2.5 at each station
d = nc.Dataset(str(CAMS_NC))
tv = nc.num2date(d['valid_time'][:], d['valid_time'].units, 'proleptic_gregorian')
YM = [f'{t.year}{t.month:02d}' for t in tv]; LAT = np.array(d['latitude'][:]); LON = np.array(d['longitude'][:]); F = np.array(d['pm2p5'][:]) * 1e9
la_ = LAT if LAT[0] < LAT[-1] else LAT[::-1]; f_ = F if LAT[0] < LAT[-1] else F[:, ::-1, :]
def bilinear(lat, lon):
    i = np.clip(np.searchsorted(la_, lat) - 1, 0, len(la_) - 2); j = np.clip(np.searchsorted(LON, lon) - 1, 0, len(LON) - 2)
    wy = (lat - la_[i]) / (la_[i + 1] - la_[i]); wx = (lon - LON[j]) / (LON[j + 1] - LON[j])
    return (1 - wy) * (1 - wx) * f_[:, i, j] + (1 - wy) * wx * f_[:, i, j + 1] + wy * (1 - wx) * f_[:, i + 1, j] + wy * wx * f_[:, i + 1, j + 1]
cams, aod, n_aod = {}, {}, 0
for sid in m.station_id.unique():
    la, lo = st.lat[sid], st.lon[sid]
    for ym, v in zip(YM, bilinear(la, lo)): cams[(sid, ym)] = v
    f = aod_file(la, lo)
    if os.path.exists(f):
        n_aod += 1
        for k, v in json.load(open(f))['properties']['parameter']['AOD_55'].items():
            if k[4:] != '13' and v is not None and v > 0: aod[(sid, k)] = v
m['cams'] = [cams.get(k) for k in zip(m.station_id, m.ym)]
m['aod'] = [aod.get(k) for k in zip(m.station_id, m.ym)]
m['l_st'] = np.log(m.cpcb); m['l_washu'] = np.log(m.sat_pm25); m['l_cams'] = np.log(m.cams); m['l_ceres'] = np.log(m.aod)
m['citykey'] = m.state + '/' + m.city.fillna(m.station_id)
print('station-months', len(m), 'stations', m.station_id.nunique(), 'with CERES', n_aod, 'cities', m.citykey.nunique())

up = m.state == 'Uttar Pradesh'
GROUPS = {'Uttar Pradesh': up, 'Rest of India': ~up, 'Rest of India without Delhi': ~up & (m.state != 'Delhi'),
          'Plain outside Uttar Pradesh': m.state.isin(IGP)}


def blocks(g, ycol, weight_city=False):
    """Per-city normal-equation blocks after demeaning within station x calendar-month cells. Base year 2019 dropped."""
    g = g.dropna(subset=[ycol]).copy()
    X = np.zeros((len(g), len(YEARS))); X[np.arange(len(g)), g.y.values - YEARS[0]] = 1.0
    cell = g.station_id + '_' + g.mo.astype(str)
    yv = g[ycol].values - g.groupby(cell)[ycol].transform('mean').values
    Xd = X - pd.DataFrame(X).groupby(cell.values).transform('mean').values
    keep = [i for i, y in enumerate(YEARS) if y != 2019]
    Xd = Xd[:, keep]
    w = (1.0 / g.groupby('citykey').station_id.transform('nunique').values) if weight_city else np.ones(len(g))
    out = {}
    for c, idx in g.groupby('citykey').indices.items():
        out[c] = ((Xd[idx] * w[idx, None]).T @ Xd[idx], (Xd[idx] * w[idx, None]).T @ yv[idx])
    return out, g.station_id.nunique(), len(g), g.groupby('y').station_id.nunique().reindex(YEARS).fillna(0).astype(int).values


def solve(A, b):
    beta = np.linalg.lstsq(A, b, rcond=None)[0]
    full = np.insert(beta, YEARS.index(2019), 0.0)                 # 2019 = 0
    return full - full[[YEARS.index(y) for y in BASE]].mean()      # relative to the 2016-2019 mean


def estimate(g, ycol, weight_city=False):
    bl, ns, nr, byyear = blocks(g, ycol, weight_city)
    keys = list(bl); A = np.array([bl[k][0] for k in keys]); B = np.array([bl[k][1] for k in keys])
    est = solve(A.sum(0), B.sum(0))
    draws = np.empty((NB, len(YEARS)))
    for i in range(NB):
        ix = rng.integers(0, len(keys), len(keys))
        draws[i] = solve(A[ix].sum(0), B[ix].sum(0))
    return est, draws, ns, nr, len(keys), byyear


rows, D = [], {}
for gname, mask in GROUPS.items():
    g = m[mask]
    for ref in ['washu', 'cams', 'ceres']:
        gg = g.dropna(subset=[f'l_{ref}']).copy(); gg['l_ratio'] = gg.l_st - gg[f'l_{ref}']
        for what, col in [('station', 'l_st'), ('reference', f'l_{ref}'), ('ratio', 'l_ratio')]:
            est, dr, ns, nr, ncity, byyear = estimate(gg, col)
            D[(gname, ref, what)] = dr
            for i, y in enumerate(YEARS):
                rows.append({'group': gname, 'reference': ref, 'series': what, 'weights': 'station-months', 'year': y, 'pct': 100 * np.expm1(est[i]),
                             'lo': 100 * np.expm1(np.percentile(dr[:, i], 2.5)), 'hi': 100 * np.expm1(np.percentile(dr[:, i], 97.5)),
                             'stations': ns, 'station_months': nr, 'cities': ncity, 'stations_in_year': int(byyear[i])})
    if gname in ('Uttar Pradesh', 'Rest of India'):
        gg = g.copy(); gg['l_ratio'] = gg.l_st - gg.l_washu
        est, dr, ns, nr, ncity, byyear = estimate(gg, 'l_ratio', weight_city=True)
        D[(gname, 'washu', 'ratio_cityw')] = dr
        for i, y in enumerate(YEARS):
            rows.append({'group': gname, 'reference': 'washu', 'series': 'ratio', 'weights': 'cities equal', 'year': y, 'pct': 100 * np.expm1(est[i]),
                         'lo': 100 * np.expm1(np.percentile(dr[:, i], 2.5)), 'hi': 100 * np.expm1(np.percentile(dr[:, i], 97.5)),
                         'stations': ns, 'station_months': nr, 'cities': ncity, 'stations_in_year': int(byyear[i])})
    print(gname, 'done', flush=True)
s = pd.DataFrame(rows); s.round(2).to_csv(f'{H}/national_series.csv', index=False)

# levels: median station / WashU ratio by year, stations with at least 6 valid months in the year (composition changes; context only)
lv = m.groupby(['station_id', 'state', 'y']).agg(n=('lr', 'size'), lr=('lr', 'mean')).reset_index()
lv = lv[lv.n >= 6]; lv['grp'] = np.where(lv.state == 'Uttar Pradesh', 'Uttar Pradesh', 'Rest of India')
first = m.groupby('station_id').y.min(); lv['entry'] = np.where(lv.station_id.map(first) <= 2019, 'reporting by 2019', 'started 2020 or later')
L1 = lv.groupby(['grp', 'y']).agg(stations=('lr', 'size'), median_ratio=('lr', lambda x: float(np.exp(np.median(x))))).reset_index()
L2 = lv.groupby(['grp', 'entry', 'y']).agg(stations=('lr', 'size'), median_ratio=('lr', lambda x: float(np.exp(np.median(x))))).reset_index()
pd.concat([L1.assign(entry='all'), L2]).round(3).to_csv(f'{H}/national_levels.csv', index=False)


def cell(gname, ref, what, years, weights='station-months'):
    x = s[(s.group == gname) & (s.reference == ref) & (s.series == what) & (s.weights == weights)].set_index('year')
    return x.loc[years]


def two_year(gname, ref, what, key=None):
    dr = D[key or (gname, ref, what)]; i = [YEARS.index(2023), YEARS.index(2024)]
    v = dr[:, i].mean(1)
    x = s[(s.group == gname) & (s.reference == ref) & (s.series == ('ratio' if what.startswith('ratio') else what)) &
          (s.weights == ('cities equal' if what == 'ratio_cityw' else 'station-months'))].set_index('year')
    est = np.log1p(x.pct[[2023, 2024]] / 100).mean()
    return 100 * np.expm1(est), 100 * np.expm1(np.percentile(v, 2.5)), 100 * np.expm1(np.percentile(v, 97.5)), v


fm = lambda v: f'{v:+.0f}%'.replace('-', '−')
L = ['# National year-by-year comparison, 2016 to 2024 (all stations)', '',
     f'{len(m):,} valid station-months, {m.station_id.nunique()} stations ({int(up.groupby(m.station_id).first().sum())} in Uttar Pradesh), {m.citykey.nunique()} cities. '
     f'Optical depth is available for {n_aod} stations. Each station is compared with itself in the same calendar month of other years. '
     'Values are relative to the 2016–2019 mean. Intervals: 95%, cities resampled.', '',
     '## Mean of 2023 and 2024 against 2016–2019', '',
     '| Group | Reference | Stations | Station PM2.5 | Reference at the stations | Station against reference |', '|---|---|---|---|---|---|']
for gname in GROUPS:
    for ref, rn in [('washu', 'WashU'), ('cams', 'CAMS'), ('ceres', 'CERES')]:
        a = two_year(gname, ref, 'station'); b = two_year(gname, ref, 'reference'); c = two_year(gname, ref, 'ratio')
        n = int(s[(s.group == gname) & (s.reference == ref) & (s.series == 'ratio')].stations.iloc[0])
        L.append(f'| {gname} | {rn} | {n} | {fm(a[0])} | {fm(b[0])} | {fm(c[0])} ({fm(c[1])} to {fm(c[2])}) |')
L += ['', '## Uttar Pradesh minus the comparison group, station against reference, mean of 2023 and 2024 (percentage points)', '',
      '| Comparison | WashU | CAMS | CERES |', '|---|---|---|---|']
for a_, b_ in [('Uttar Pradesh', 'Rest of India'), ('Uttar Pradesh', 'Rest of India without Delhi'), ('Uttar Pradesh', 'Plain outside Uttar Pradesh')]:
    cells = []
    for ref in ['washu', 'cams', 'ceres']:
        x = two_year(a_, ref, 'ratio'); y = two_year(b_, ref, 'ratio')
        dd = 100 * np.expm1(x[3]) - 100 * np.expm1(y[3])
        cells.append(f'{x[0] - y[0]:+.0f} ({np.percentile(dd, 2.5):+.0f} to {np.percentile(dd, 97.5):+.0f})'.replace('-', '−'))
    L.append(f'| {a_} minus {b_} | ' + ' | '.join(cells) + ' |')
xw = two_year('Uttar Pradesh', 'washu', 'ratio_cityw', ('Uttar Pradesh', 'washu', 'ratio_cityw')); yw = two_year('Rest of India', 'washu', 'ratio_cityw', ('Rest of India', 'washu', 'ratio_cityw'))
L += ['', f'Cities weighted equally, against WashU: Uttar Pradesh {fm(xw[0])} ({fm(xw[1])} to {fm(xw[2])}); rest of India {fm(yw[0])} ({fm(yw[1])} to {fm(yw[2])}).', '',
      '## By year, against WashU', '', '| Year | UP: stations in year | UP station PM2.5 | UP WashU | UP station against WashU | Rest: stations in year | Rest station PM2.5 | Rest WashU | Rest station against WashU |', '|---|---|---|---|---|---|---|---|---|']
for y in YEARS:
    c = lambda g, w: cell(g, 'washu', w, y)
    u1, u2, u3 = c('Uttar Pradesh', 'station'), c('Uttar Pradesh', 'reference'), c('Uttar Pradesh', 'ratio')
    r1, r2, r3 = c('Rest of India', 'station'), c('Rest of India', 'reference'), c('Rest of India', 'ratio')
    L.append(f'| {y} | {int(u3.stations_in_year)} | {fm(u1.pct)} | {fm(u2.pct)} | {fm(u3.pct)} ({fm(u3.lo)} to {fm(u3.hi)}) | {int(r3.stations_in_year)} | {fm(r1.pct)} | {fm(r2.pct)} | {fm(r3.pct)} ({fm(r3.lo)} to {fm(r3.hi)}) |')
L += ['', '## Level: median ratio of station PM2.5 to WashU, by year (stations with at least six valid months; the set of stations changes)', '',
      '| Year | UP, all (n) | UP, reporting by 2019 (n) | UP, started 2020 or later (n) | Rest, all (n) | Rest, reporting by 2019 (n) | Rest, started 2020 or later (n) |', '|---|---|---|---|---|---|---|']
def lvl(grp, entry, y):
    t = (L1 if entry == 'all' else L2[L2.entry == entry]); t = t[(t.grp == grp) & (t.y == y)]
    return '' if t.empty else f'{t.median_ratio.iloc[0]:.2f} ({int(t.stations.iloc[0])})'
for y in YEARS:
    L.append(f'| {y} | ' + ' | '.join(lvl(g, e, y) for g in ['Uttar Pradesh', 'Rest of India'] for e in ['all', 'reporting by 2019', 'started 2020 or later']) + ' |')
open(f'{H}/national_series.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
