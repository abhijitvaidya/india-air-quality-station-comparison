#!/usr/bin/env python3
"""PM10 against CAMS PM10: checks on the base period and on the manual stations (Section 4.2, 4.6; Table 4; Appendix D1).
A. Continuous PM10, Uttar Pradesh minus rest of India, 2023-2024: interval with all draws, with only the draws in which
   every base year has data in both groups, and with later bases (2017-2019, 2018-2019, 2019).
B. Manual stations from 2019: all cities; cities that also have continuous stations; cities with manual stations only;
   with the station-months flagged as possible reprints kept.
C. By year, manual and continuous stations in the cities with both networks.

Run:  python code/20_pm10_manual_checks.py
Output: results/pm10_manual_checks.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, numpy as np, pandas as pd, netCDF4 as nc
H = str(RESULTS); DS = str(DERIVED); NB = 2000
rng = np.random.default_rng(20261124)
d = nc.Dataset(str(CAMS_NC))
tv = nc.num2date(d['valid_time'][:], d['valid_time'].units, 'proleptic_gregorian')
YM = [f'{t.year}{t.month:02d}' for t in tv]; LAT = np.array(d['latitude'][:]); LON = np.array(d['longitude'][:]); F = np.array(d['pm10'][:]) * 1e9
la_ = LAT if LAT[0] < LAT[-1] else LAT[::-1]; f_ = F if LAT[0] < LAT[-1] else F[:, ::-1, :]
def bil(lat, lon):
    i = np.clip(np.searchsorted(la_, lat) - 1, 0, len(la_) - 2); j = np.clip(np.searchsorted(LON, lon) - 1, 0, len(LON) - 2)
    wy = (lat - la_[i]) / (la_[i + 1] - la_[i]); wx = (lon - LON[j]) / (LON[j + 1] - LON[j])
    return (1 - wy) * (1 - wx) * f_[:, i, j] + (1 - wy) * wx * f_[:, i, j + 1] + wy * (1 - wx) * f_[:, i + 1, j] + wy * wx * f_[:, i + 1, j + 1]
st = pd.read_csv(f'{DS}/stations.csv').set_index('station_id'); st = st[st.lat.notna()]
c = pd.read_csv(f'{DS}/continuous_station_month_pm10.csv', dtype={'ym': str}); c = c[(c.ndays >= 15) & (c['mean'] > 0) & c.station_id.isin(st.index)].copy()
c['y'] = c.ym.str[:4].astype(int); c['mo'] = c.ym.str[4:].astype(int); c = c[(c.y >= 2016) & (c.y <= 2025)]
c['state'] = c.station_id.map(st.state); c['city'] = c.station_id.map(st.city); c['citykey'] = c.state + '/' + c.city.fillna(c.station_id)
cams = {}
for sid in c.station_id.unique():
    for ym, v in zip(YM, bil(st.lat[sid], st.lon[sid])): cams[(sid, ym)] = v
c['lr'] = np.log(c['mean']) - np.log([cams[k] for k in zip(c.station_id, c.ym)])
up = c.state == 'Uttar Pradesh'
CALIAS = {'Prayagraj': ['Prayagraj', 'Allahabad'], 'Raebareli': ['Raebareli', 'Rae Bareli']}
geo = pd.read_csv(f'{DS}/manual_city_centres.csv').set_index('city')
def manual(keep_reprints):
    m = pd.read_csv(f'{DS}/manual_station_month.csv'); m = m[m.pm10.notna() & (m.pm10 > 0)]
    if not keep_reprints: m = m[m.repeat_of.isna()]
    m = m.copy(); m['y'] = m.year.astype(int); m['mo'] = m.month.astype(int); m['ym'] = m.y.astype(str) + m.mo.map('{:02d}'.format)
    m['station_id'] = 'manual/' + m.city + '/' + m.station; m['citykey'] = 'UP/' + m.city
    xy, both = {}, set()
    for city in m.city.unique():
        ids = [s for s in st.index if st.state[s] == 'Uttar Pradesh' and st.city[s] in CALIAS.get(city, [city])]
        if ids: xy[city] = (st.lat[ids].mean(), st.lon[ids].mean()); both.add(city)
        else: xy[city] = (geo.lat[city], geo.lon[city])
    cc = {city: dict(zip(YM, bil(*p))) for city, p in xy.items()}
    m['lr'] = np.log(m.pm10) - np.log([cc[ci][ym] for ci, ym in zip(m.city, m.ym)]); return m, both


def estimate(g, years, base):
    g = g[g.y.isin(years)]; ref = base[-1]; keep = [i for i, y in enumerate(years) if y != ref]; bi = [years.index(y) for y in base]
    X = np.zeros((len(g), len(years))); X[np.arange(len(g)), g.y.values - years[0]] = 1.0
    cell = (g.station_id + '_' + g.mo.astype(str)).values
    yv = g.lr.values - g.groupby(cell).lr.transform('mean').values
    Xd = (X - pd.DataFrame(X).groupby(cell).transform('mean').values)[:, keep]
    keys, inv = np.unique(g.citykey.values, return_inverse=True)
    A = np.zeros((len(keys), len(keep), len(keep))); B = np.zeros((len(keys), len(keep)))
    for k in range(len(keys)):
        idx = inv == k; A[k] = Xd[idx].T @ Xd[idx]; B[k] = Xd[idx].T @ yv[idx]
    nref = np.array([((inv == k) & (g.y.values == ref)).sum() for k in range(len(keys))], float)
    def solve(As, Bs):
        dead = np.einsum('...ii->...i', As) < 1e-9
        beta = np.where(dead, np.nan, np.einsum('...ij,...j->...i', np.linalg.pinv(As, hermitian=True), Bs))
        full = np.insert(beta, years.index(ref), 0.0, axis=-1); miss = np.isnan(full[..., bi]).any(-1)
        return full - np.nanmean(full[..., bi], axis=-1, keepdims=True), miss
    est, _ = solve(A.sum(0), B.sum(0))
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    dr, miss = solve(np.tensordot(W, A, 1), W @ B); bad = (W @ nref) == 0; dr[bad] = np.nan
    return est, dr, miss | bad, g.station_id.nunique(), len(keys)


pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
L = ['# PM10 against CAMS PM10: base period and manual stations', '', '## A. Continuous PM10 against CAMS PM10: Uttar Pradesh minus rest of India, mean of 2023 and 2024', '',
     '| Base years | Draws used | UP | Rest | Difference (points) | 95% interval | Draws with a base year missing in UP |', '|---|---|---|---|---|---|---|']
Y = list(range(2016, 2025))
for base, yrs in [([2016, 2017, 2018, 2019], Y), ([2017, 2018, 2019], Y[1:]), ([2018, 2019], Y[2:]), ([2019], Y[3:])]:
    eu, du, mu, _, _ = estimate(c[up], yrs, base); er, dr_, mr, _, _ = estimate(c[~up], yrs, base)
    i = [yrs.index(2023), yrs.index(2024)]; a, b = pc(eu[i].mean()), pc(er[i].mean()); dd = pc(du[:, i].mean(1)) - pc(dr_[:, i].mean(1))
    for label, ok in [('all', np.ones(NB, bool)), ('complete base only', ~(mu | mr))]:
        lo, hi = np.nanpercentile(dd[ok], [2.5, 97.5])
        L.append(f'| {base[0]}–{base[-1]} | {label} ({int(ok.sum())}) | {fm(a)}% | {fm(b)}% | {fm(a - b)} | {fm(lo)} to {fm(hi)} | {100 * mu.mean():.0f}% |')
L += ['', '## B. Manual and continuous PM10 against CAMS PM10, relative to 2019', '',
      '| Group | Stations (cities) | 2022 | 2024 | 2023–24 (interval) | 2025 |', '|---|---|---|---|---|---|']
YB = list(range(2019, 2026)); m0, both = manual(False); m1, _ = manual(True)
inM = up & c.city.isin(sum((CALIAS.get(x, [x]) for x in m0.city.unique()), []))
for name, g in [('Manual, all cities', m0), ('Manual, cities with both networks', m0[m0.city.isin(both)]), ('Manual, cities with manual stations only', m0[~m0.city.isin(both)]),
                ('Continuous, cities with both networks', c[inM]), ('Manual, all cities, flagged reprints kept', m1), ('Manual, cities with both networks, flagged reprints kept', m1[m1.city.isin(both)])]:
    e, dr, _, ns, nk = estimate(g, YB, [2019]); i = [YB.index(2023), YB.index(2024)]; v = dr[:, i].mean(1)
    L.append(f'| {name} | {ns} ({nk}) | {fm(pc(e[YB.index(2022)]))}% | {fm(pc(e[YB.index(2024)]))}% | {fm(pc(e[i].mean()))}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}) | {fm(pc(e[YB.index(2025)]))}% |')
L += ['', f'Cities with both networks: {len(both)}. Manual stations reporting in 2024: {m0[m0.y == 2024].station_id.nunique()} in {m0[m0.y == 2024].city.nunique()} cities; in cities with both networks: {m0[(m0.y == 2024) & m0.city.isin(both)].station_id.nunique()}.']
L += ['', '## C. By year (stations), relative to 2019', '', '| Year | Manual, cities with both networks | Continuous, cities with both networks |', '|---|---|---|']
ea, _, _, _, _ = estimate(m0[m0.city.isin(both)], YB, [2019]); eb, _, _, _, _ = estimate(c[inM], YB, [2019])
na = m0[m0.city.isin(both)].groupby('y').station_id.nunique(); nb = c[inM].groupby('y').station_id.nunique()
for i, y in enumerate(YB): L.append(f'| {y} | {fm(pc(ea[i]))}% ({na.get(y, 0)}) | {fm(pc(eb[i]))}% ({nb.get(y, 0)}) |')
# PM10 on the 2018-2019 base: group values with intervals
yrs = Y[2:]; i = [yrs.index(2023), yrs.index(2024)]
for nm, g in [('Uttar Pradesh', c[up]), ('Rest of India', c[~up])]:
    e, dr, _, ns, nk = estimate(g, yrs, [2018, 2019]); v = dr[:, i].mean(1)
    L.append(f'PM10 against CAMS, base 2018–2019, {nm}: {fm(pc(e[i].mean()))}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}); stations {ns}; UP stations by year ' + str(g.groupby("y").station_id.nunique().to_dict()) if nm == 'Uttar Pradesh' else f'PM10 against CAMS, base 2018–2019, {nm}: {fm(pc(e[i].mean()))}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}); stations {ns}')
open(f'{H}/pm10_manual_checks.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
