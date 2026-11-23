#!/usr/bin/env python3
"""PM10 against CAMS PM10: continuous stations nationally, and manual and continuous stations in Uttar Pradesh
(Section 3.5, Table 4).

Method as 11_national_series.py: each station x calendar-month cell has its own mean, year effects are estimated by least
squares on cell-demeaned logs, and the same estimator is applied to log(station), log(reference) and log(station /
reference) on identical station-months. Intervals: cluster bootstrap, 2,000 draws, percentile 95%. A year with no data in
a draw is left out of that draw.

Part A. Continuous PM10, 2016-2024, relative to 2016-2019 (Table A3 uses 21_pm10_base_periods.py for other bases).
Part B. From 2019 to 2025, relative to 2019, against CAMS PM10 only. Manual stations: all stations in the UPPCB monthly
  reports, station-months flagged as possible reprints left out. A manual station uses CAMS at the mean location of its
  city's continuous stations, or at the city centre where the city has none. Clusters: cities for the groups, stations
  inside a single city.

Run:  python code/19_pm10_manual.py
Outputs: results/pm10_manual_continuous.csv, results/pm10_manual_continuous.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os
import numpy as np, pandas as pd
import netCDF4 as nc

H = str(RESULTS); DS = str(DERIVED); NB = 2000
rng = np.random.default_rng(20261121)
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']

# ---------- CAMS PM10 ----------
d = nc.Dataset(str(CAMS_NC))
tv = nc.num2date(d['valid_time'][:], d['valid_time'].units, 'proleptic_gregorian')
YM = [f'{t.year}{t.month:02d}' for t in tv]; LAT = np.array(d['latitude'][:]); LON = np.array(d['longitude'][:]); F = np.array(d['pm10'][:]) * 1e9
la_ = LAT if LAT[0] < LAT[-1] else LAT[::-1]; f_ = F if LAT[0] < LAT[-1] else F[:, ::-1, :]
def bilinear(lat, lon):
    i = np.clip(np.searchsorted(la_, lat) - 1, 0, len(la_) - 2); j = np.clip(np.searchsorted(LON, lon) - 1, 0, len(LON) - 2)
    wy = (lat - la_[i]) / (la_[i + 1] - la_[i]); wx = (lon - LON[j]) / (LON[j + 1] - LON[j])
    return (1 - wy) * (1 - wx) * f_[:, i, j] + (1 - wy) * wx * f_[:, i, j + 1] + wy * (1 - wx) * f_[:, i + 1, j] + wy * wx * f_[:, i + 1, j + 1]

# ---------- continuous PM10 station-months ----------
st = pd.read_csv(f'{DS}/stations.csv').set_index('station_id')
c = pd.read_csv(f'{DS}/continuous_station_month_pm10.csv', dtype={'ym': str})
c = c[(c.ndays >= 15) & (c['mean'] > 0)].rename(columns={'mean': 'pm10'})
c = c[c.station_id.isin(st.index)].copy()
c['y'] = c.ym.str[:4].astype(int); c['mo'] = c.ym.str[4:].astype(int)
c = c[(c.y >= 2016) & (c.y <= 2025)]
c['state'] = c.station_id.map(st.state); c['city'] = c.station_id.map(st.city)
c['citykey'] = c.state + '/' + c.city.fillna(c.station_id)
cams, aod, n_aod = {}, {}, 0
for sid in c.station_id.unique():
    la, lo = st.lat[sid], st.lon[sid]
    for ym, v in zip(YM, bilinear(la, lo)): cams[(sid, ym)] = v
    f = aod_file(la, lo)
    if os.path.exists(f):
        n_aod += 1
        for k, v in json.load(open(f))['properties']['parameter']['AOD_55'].items():
            if k[4:] != '13' and v is not None and v > 0: aod[(sid, k)] = v
sat = pd.read_csv(f'{DS}/satellite_pm25_station_month_v6.csv', dtype={'ym': str}).dropna()
sat = sat[sat.sat_pm25 > 0].set_index(['station_id', 'ym']).sat_pm25
c['cams'] = [cams.get(k) for k in zip(c.station_id, c.ym)]
c['aod'] = [aod.get(k) for k in zip(c.station_id, c.ym)]
c['washu'] = [sat.get(k, np.nan) for k in zip(c.station_id, c.ym)]
c['l_st'] = np.log(c.pm10)
for r_ in ['cams', 'aod', 'washu']: c[f'l_{r_}'] = np.log(c[r_].astype(float))
c = c.rename(columns={'l_aod': 'l_ceres'})
print('continuous PM10 station-months', len(c), 'stations', c.station_id.nunique(), 'with CERES cache', n_aod, 'cities', c.citykey.nunique(), flush=True)

# ---------- manual PM10 station-months ----------
man = pd.read_csv(f'{DS}/manual_station_month.csv')
man = man[man.pm10.notna() & (man.pm10 > 0)]; man = man[man.repeat_of.isna()].copy()
man['y'] = man.year.astype(int); man['mo'] = man.month.astype(int); man['ym'] = man.y.astype(str) + man.mo.map('{:02d}'.format)
man['station_id'] = 'manual/' + man.city + '/' + man.station; man['citykey'] = 'Uttar Pradesh/' + man.city; man['state'] = 'Uttar Pradesh'
CALIAS = {'Prayagraj': ['Prayagraj', 'Allahabad'], 'Raebareli': ['Raebareli', 'Rae Bareli']}
geo = pd.read_csv(f'{DS}/manual_city_centres.csv').set_index('city')
cityxy, cityxy_src = {}, {}
for city in sorted(man.city.unique()):
    ids = [s for s in st.index if st.state[s] == 'Uttar Pradesh' and st.city[s] in CALIAS.get(city, [city]) and pd.notna(st.lat[s])]
    if ids: cityxy[city] = (st.lat[ids].mean(), st.lon[ids].mean()); cityxy_src[city] = f'{len(ids)} continuous stations'
    else: cityxy[city] = (geo.lat[city], geo.lon[city]); cityxy_src[city] = 'city centre (OSM)'
ccity = {city: dict(zip(YM, bilinear(*xy))) for city, xy in cityxy.items()}
man['cams'] = [ccity[ci].get(ym) for ci, ym in zip(man.city, man.ym)]
man = man.dropna(subset=['cams']); man['l_st'] = np.log(man.pm10); man['l_cams'] = np.log(man.cams)
MCITIES = sorted(man.city.unique())
print('manual PM10 station-months', len(man), 'stations', man.station_id.nunique(), 'cities', len(MCITIES), 'years', man.y.min(), man.y.max(), flush=True)


def estimate(g, ycol, years, base, cluster):
    """Year effects relative to the mean of the base years; cluster bootstrap. Returns est, draws, n stations, n rows, n clusters, stations by year."""
    g = g.dropna(subset=[ycol]); g = g[g.y.isin(years)]
    y0 = years[0]; ref = base[-1]; keep = [i for i, y in enumerate(years) if y != ref]
    X = np.zeros((len(g), len(years))); X[np.arange(len(g)), g.y.values - y0] = 1.0
    cell = (g.station_id + '_' + g.mo.astype(str)).values
    yv = g[ycol].values - g.groupby(cell)[ycol].transform('mean').values
    Xd = (X - pd.DataFrame(X).groupby(cell).transform('mean').values)[:, keep]
    keys, inv = np.unique(g[cluster].values, return_inverse=True)
    A = np.zeros((len(keys), len(keep), len(keep))); B = np.zeros((len(keys), len(keep)))
    for k in range(len(keys)):
        idx = inv == k; A[k] = Xd[idx].T @ Xd[idx]; B[k] = Xd[idx].T @ yv[idx]
    bi = [years.index(y) for y in base]
    nref = np.array([((inv == k) & (g.y.values == ref)).sum() for k in range(len(keys))], float)   # rows in the reference year

    def solve(As, Bs):
        dead = np.einsum('...ii->...i', As) < 1e-9                      # year without data in this sample
        beta = np.einsum('...ij,...j->...i', np.linalg.pinv(As, hermitian=True), Bs)
        beta = np.where(dead, np.nan, beta)
        full = np.insert(beta, years.index(ref), 0.0, axis=-1)
        return full - np.nanmean(full[..., bi], axis=-1, keepdims=True)
    est = solve(A.sum(0), B.sum(0))
    if nref.sum() == 0: est = np.full(len(years), np.nan)            # no data in the reference year: not identified
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    draws = solve(np.tensordot(W, A, 1), W @ B)
    draws[(W @ nref) == 0] = np.nan
    byyear = g.groupby('y').station_id.nunique().reindex(years).fillna(0).astype(int).values
    return est, draws, g.station_id.nunique(), len(g), len(keys), byyear


rows, D = [], {}
def run(part, gname, g, refs, years, base, cluster='citykey'):
    for ref in refs:
        gg = g.dropna(subset=[f'l_{ref}']).copy(); gg['l_ratio'] = gg.l_st - gg[f'l_{ref}']
        if gg.empty: continue
        for what, col in [('station', 'l_st'), ('reference', f'l_{ref}'), ('ratio', 'l_ratio')]:
            est, dr, ns, nr, ncl, byyear = estimate(gg, col, years, base, cluster)
            D[(part, gname, ref, what)] = (est, dr, years)
            for i, y in enumerate(years):
                ok = np.isfinite(dr[:, i])
                rows.append({'part': part, 'group': gname, 'reference': ref, 'series': what, 'year': y, 'pct': 100 * np.expm1(est[i]),
                             'lo': 100 * np.expm1(np.nanpercentile(dr[:, i], 2.5)) if ok.any() else np.nan,
                             'hi': 100 * np.expm1(np.nanpercentile(dr[:, i], 97.5)) if ok.any() else np.nan,
                             'stations': ns, 'station_months': nr, 'clusters': ncl, 'cluster_unit': cluster, 'stations_in_year': int(byyear[i])})
    print(part, gname, 'done', flush=True)


# ---------- Part A ----------
YA = list(range(2016, 2025)); BA = [2016, 2017, 2018, 2019]
up = c.state == 'Uttar Pradesh'; cA = c[c.y <= 2024]
GA = {'Uttar Pradesh': up, 'Rest of India': ~up, 'Rest of India without Delhi': ~up & (c.state != 'Delhi'), 'Plain outside Uttar Pradesh': c.state.isin(IGP)}
for gname, mask in GA.items(): run('A', gname, cA[mask.reindex(cA.index)], ['cams', 'washu', 'ceres'], YA, BA)

# ---------- Part B ----------
YB = list(range(2019, 2026)); BB = [2019]
inM = up & c.city.isin(sum((CALIAS.get(x, [x]) for x in MCITIES), []))
run('B', 'Manual stations, all Uttar Pradesh reports', man, ['cams'], YB, BB)
run('B', 'Continuous stations, cities with manual stations', c[inM], ['cams'], YB, BB)
run('B', 'Manual stations, cities without continuous stations', man[man.city.map(lambda x: cityxy_src[x].startswith('city'))], ['cams'], YB, BB)
run('B', 'Continuous stations, Uttar Pradesh', c[up], ['cams'], YB, BB)
run('B', 'Continuous stations, rest of India', c[~up], ['cams'], YB, BB)
for city in ['Varanasi']:
    run('B', f'{city}: manual', man[man.city == city], ['cams'], YB, BB, cluster='station_id')
    run('B', f'{city}: continuous', c[up & (c.city == city)], ['cams'], YB, BB, cluster='station_id')

s = pd.DataFrame(rows); s.round(2).to_csv(f'{H}/pm10_manual_continuous.csv', index=False)


def span(part, gname, ref, what, yrs):
    """Mean of the given years (log scale) with interval."""
    if (part, gname, ref, what) not in D: return None
    est, dr, years = D[(part, gname, ref, what)]; i = [years.index(y) for y in yrs]
    v = dr[:, i].mean(1); e = est[i].mean()
    if not np.isfinite(v).any(): return np.nan, np.nan, np.nan, v
    return 100 * np.expm1(e), 100 * np.expm1(np.nanpercentile(v, 2.5)), 100 * np.expm1(np.nanpercentile(v, 97.5)), v


fm = lambda v: 'n/a' if v is None or not np.isfinite(v) else f'{v:+.0f}%'.replace('-', '−')
def n_of(part, gname, ref):
    x = s[(s.part == part) & (s.group == gname) & (s.reference == ref) & (s.series == 'ratio')]
    return int(x.stations.iloc[0]) if len(x) else 0
RN = {'cams': 'CAMS PM10', 'washu': 'WashU PM2.5', 'ceres': 'CERES AOD'}
L = ['# PM10 against CAMS PM10, by year', '',
     f'Continuous PM10: {len(cA):,} valid station-months in 2016–2024 at {cA.station_id.nunique()} stations; the {cA.dropna(subset=["l_cams"]).station_id.nunique()} with coordinates are compared. '
     f'Manual PM10: {len(man):,} monthly values, {man.station_id.nunique()} stations in {len(MCITIES)} Uttar Pradesh cities, {man.y.min()}–{man.y.max()} (CAMS at the city centre for {sum(v.startswith("city") for v in cityxy_src.values())} cities without continuous stations). '
     'Each station is compared with itself in the same calendar month of other years. Intervals: 95%.', '',
     '## A. Continuous PM10, mean of 2023 and 2024 against 2016–2019', '',
     '| Group | Reference | Stations | Station PM10 | Reference at the stations | Station against reference |', '|---|---|---|---|---|---|']
for gname in GA:
    for ref in ['cams', 'washu', 'ceres']:
        a, b, r = (span('A', gname, ref, w, [2023, 2024]) for w in ['station', 'reference', 'ratio'])
        if r is None: continue
        L.append(f'| {gname} | {RN[ref]} | {n_of("A", gname, ref)} | {fm(a[0])} | {fm(b[0])} | {fm(r[0])} ({fm(r[1])} to {fm(r[2])}) |')
L += ['', '## A. Uttar Pradesh minus the comparison group, station against reference (percentage points)', '', '| Comparison | CAMS PM10 | WashU PM2.5 | CERES AOD |', '|---|---|---|---|']
for b_ in ['Rest of India', 'Rest of India without Delhi', 'Plain outside Uttar Pradesh']:
    cells = []
    for ref in ['cams', 'washu', 'ceres']:
        x = span('A', 'Uttar Pradesh', ref, 'ratio', [2023, 2024]); y = span('A', b_, ref, 'ratio', [2023, 2024])
        dd = 100 * np.expm1(x[3]) - 100 * np.expm1(y[3])
        cells.append(f'{x[0] - y[0]:+.0f} ({np.nanpercentile(dd, 2.5):+.0f} to {np.nanpercentile(dd, 97.5):+.0f})'.replace('-', '−'))
    L.append(f'| Uttar Pradesh minus {b_} | ' + ' | '.join(cells) + ' |')
L += ['', '## A. By year, continuous PM10 against CAMS PM10', '',
      '| Year | UP: stations in year | UP station PM10 | UP CAMS PM10 | UP station against CAMS | Rest: stations in year | Rest station PM10 | Rest CAMS PM10 | Rest station against CAMS |', '|---|---|---|---|---|---|---|---|---|']
def cellv(part, g, w, y, ref='cams'):
    return s[(s.part == part) & (s.group == g) & (s.reference == ref) & (s.series == w) & (s.year == y)].iloc[0]
for y in YA:
    u1, u2, u3 = (cellv('A', 'Uttar Pradesh', w, y) for w in ['station', 'reference', 'ratio'])
    r1, r2, r3 = (cellv('A', 'Rest of India', w, y) for w in ['station', 'reference', 'ratio'])
    L.append(f'| {y} | {int(u3.stations_in_year)} | {fm(u1.pct)} | {fm(u2.pct)} | {fm(u3.pct)} ({fm(u3.lo)} to {fm(u3.hi)}) | {int(r3.stations_in_year)} | {fm(r1.pct)} | {fm(r2.pct)} | {fm(r3.pct)} ({fm(r3.lo)} to {fm(r3.hi)}) |')

L += ['', '## B. From 2019: manual and continuous PM10 against CAMS PM10, relative to 2019', '',
      'Manual stations use CAMS PM10 at the mean location of the city\'s continuous stations, or at the city centre.', '',
      '| Group | Stations | Station PM10, 2023–24 | CAMS PM10, 2023–24 | Station against CAMS, 2023–24 | Station against CAMS, 2025 |', '|---|---|---|---|---|---|']
GB = ['Manual stations, all Uttar Pradesh reports', 'Continuous stations, cities with manual stations', 'Continuous stations, rest of India',
      'Manual stations, cities without continuous stations',
      'Continuous stations, Uttar Pradesh', 'Varanasi: manual', 'Varanasi: continuous']
for gname in GB:
    a, b, r, r25 = (span('B', gname, 'cams', w, yy) for w, yy in [('station', [2023, 2024]), ('reference', [2023, 2024]), ('ratio', [2023, 2024]), ('ratio', [2025])])
    if r is None: L.append(f'| {gname} | 0 | n/a | n/a | n/a | n/a |'); continue
    L.append(f'| {gname} | {n_of("B", gname, "cams")} | {fm(a[0])} | {fm(b[0])} | {fm(r[0])} ({fm(r[1])} to {fm(r[2])}) | {fm(r25[0])} ({fm(r25[1])} to {fm(r25[2])}) |')
L += ['', '## B. By year, station against CAMS PM10, relative to 2019', '', '| Year | Manual, all reports (n) | Continuous, cities with manual stations (n) | Continuous, rest of India (n) |', '|---|---|---|---|']
for y in YB:
    cs = []
    for gname in GB[:3]:
        x = cellv('B', gname, 'ratio', y); cs.append(f'{fm(x.pct)} ({int(x.stations_in_year)})')
    L.append(f'| {y} | ' + ' | '.join(cs) + ' |')
open(f'{H}/pm10_manual_continuous.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
