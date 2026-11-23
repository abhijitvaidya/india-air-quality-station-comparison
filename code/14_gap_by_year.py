#!/usr/bin/env python3
"""Uttar Pradesh minus the rest of India, year by year, against each reference, with city-bootstrap intervals (Table E2).
Also for the 166 eligible stations of the single-station comparison, and leaving out one year at a time.

Run after 13_single_station.py:  python code/14_gap_by_year.py
Outputs: results/gap_by_year.csv, results/gap_by_year.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys
import numpy as np, pandas as pd, netCDF4 as nc
H = str(RESULTS)
from common import paired
import two_period as E
YEARS = list(range(2016, 2025)); BASE = [2016, 2017, 2018, 2019]; REF = 2019; NB = 2000
rng = np.random.default_rng(20261123)
m = paired('v6'); st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id')
d = nc.Dataset(str(CAMS_NC))
tv = nc.num2date(d['valid_time'][:], d['valid_time'].units, 'proleptic_gregorian')
YM = [f'{t.year}{t.month:02d}' for t in tv]; LAT = np.array(d['latitude'][:]); LON = np.array(d['longitude'][:]); F = np.array(d['pm2p5'][:]) * 1e9
la_ = LAT if LAT[0] < LAT[-1] else LAT[::-1]; f_ = F if LAT[0] < LAT[-1] else F[:, ::-1, :]
def bil(lat, lon):
    i = np.clip(np.searchsorted(la_, lat) - 1, 0, len(la_) - 2); j = np.clip(np.searchsorted(LON, lon) - 1, 0, len(LON) - 2)
    wy = (lat - la_[i]) / (la_[i + 1] - la_[i]); wx = (lon - LON[j]) / (LON[j + 1] - LON[j])
    return (1 - wy) * (1 - wx) * f_[:, i, j] + (1 - wy) * wx * f_[:, i, j + 1] + wy * (1 - wx) * f_[:, i + 1, j] + wy * wx * f_[:, i + 1, j + 1]
cams, aod = {}, {}
for sid in m.station_id.unique():
    la, lo = st.lat[sid], st.lon[sid]
    for ym, v in zip(YM, bil(la, lo)): cams[(sid, ym)] = v
    f = aod_file(la, lo)
    if os.path.exists(f):
        for k, v in json.load(open(f))['properties']['parameter']['AOD_55'].items():
            if k[4:] != '13' and v is not None and v > 0: aod[(sid, k)] = v
m['l_washu'] = np.log(m.sat_pm25)
m['l_cams'] = np.log(pd.Series([cams.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float))
m['l_ceres'] = np.log(pd.Series([aod.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float))
m['l_st'] = np.log(m.cpcb); m['citykey'] = m.state + '/' + m.city.fillna(m.station_id); up = m.state == 'Uttar Pradesh'
bi = [YEARS.index(y) for y in BASE]; keep = [i for i, y in enumerate(YEARS) if y != REF]


def estimate(g, ycol):
    g = g.dropna(subset=[ycol]); X = np.zeros((len(g), len(YEARS))); X[np.arange(len(g)), g.y.values - YEARS[0]] = 1.0
    cell = (g.station_id + '_' + g.mo.astype(str)).values
    yv = g[ycol].values - g.groupby(cell)[ycol].transform('mean').values
    Xd = (X - pd.DataFrame(X).groupby(cell).transform('mean').values)[:, keep]
    keys, inv = np.unique(g.citykey.values, return_inverse=True)
    A = np.zeros((len(keys), len(keep), len(keep))); B = np.zeros((len(keys), len(keep)))
    for k in range(len(keys)):
        idx = inv == k; A[k] = Xd[idx].T @ Xd[idx]; B[k] = Xd[idx].T @ yv[idx]
    nref = np.array([((inv == k) & (g.y.values == REF)).sum() for k in range(len(keys))], float)
    def solve(As, Bs):
        dead = np.einsum('...ii->...i', As) < 1e-9
        beta = np.where(dead, np.nan, np.einsum('...ij,...j->...i', np.linalg.pinv(As, hermitian=True), Bs))
        full = np.insert(beta, YEARS.index(REF), 0.0, axis=-1)
        return full - np.nanmean(full[..., bi], axis=-1, keepdims=True)
    est = solve(A.sum(0), B.sum(0))
    W = rng.multinomial(len(keys), np.full(len(keys), 1 / len(keys)), size=NB).astype(float)
    dr = solve(np.tensordot(W, A, 1), W @ B); dr[(W @ nref) == 0] = np.nan
    return est, dr, g.station_id.nunique()


pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
RN = {'washu': 'WashU', 'cams': 'CAMS', 'ceres': 'CERES'}
panel = set(pd.read_csv(f'{H}/single_station.csv').station_id)
rows = []; L = ['# Uttar Pradesh minus rest of India, by year, with intervals', '',
                'Station-to-reference ratio, relative to the 2016–2019 mean; difference in percentage points; 95% intervals, cities resampled.', '']
for label, mask in [('All stations', pd.Series(True, index=m.index)), ('166 stations of the single-station comparison', m.station_id.isin(panel))]:
    L += [f'## {label}', '', '| Year | ' + ' | '.join(f'{RN[r]}: UP, rest, difference' for r in RN) + ' |', '|---|---|---|---|']
    res = {}
    for ref in RN:
        mm = m[mask].dropna(subset=[f'l_{ref}']).copy(); mm['lr_'] = mm.l_st - mm[f'l_{ref}']
        res[ref] = (estimate(mm[mm.state == 'Uttar Pradesh'], 'lr_'), estimate(mm[mm.state != 'Uttar Pradesh'], 'lr_'))
    def cellf(ref, idx):
        (eu, du, nu), (er, dr_, nr) = res[ref]
        a = pc(eu[idx].mean()); b = pc(er[idx].mean()); dd = pc(du[:, idx].mean(1)) - pc(dr_[:, idx].mean(1))
        lo, hi = np.nanpercentile(dd, [2.5, 97.5]); return a, b, a - b, lo, hi, nu, nr
    for y in YEARS + ['2023–2024']:
        idx = [YEARS.index(2023), YEARS.index(2024)] if y == '2023–2024' else [YEARS.index(y)]
        cs = []
        for ref in RN:
            a, b, dff, lo, hi, nu, nr = cellf(ref, idx)
            cs.append(f'{fm(a)}%, {fm(b)}%, {fm(dff)} ({fm(lo)} to {fm(hi)})')
            rows.append({'stations': label, 'reference': ref, 'year': y, 'up_pct': a, 'rest_pct': b, 'diff_points': dff, 'lo': lo, 'hi': hi, 'n_up': nu, 'n_rest': nr})
        L.append(f'| {y} | ' + ' | '.join(cs) + ' |')
    L += ['', f'Stations: Uttar Pradesh {res["washu"][0][2]}, rest {res["washu"][1][2]}.', '']
pd.DataFrame(rows).round(2).to_csv(f'{H}/gap_by_year.csv', index=False)

# leave-one-year-out on the 166 stations, against WashU
L += ['## Leave one year out, 166 stations, against WashU', '', '| Year left out | Stations with a result | Ardhali Bazar change | Rank | Next station |', '|---|---|---|---|---|']
mp = m[m.station_id.isin(panel)].copy(); mp['lr'] = mp.l_st - mp.l_washu
for yr in [None, 2016, 2017, 2018, 2019, 2023, 2024]:
    out = {}
    for sid, g in mp[mp.y != yr].groupby('station_id') if yr else mp.groupby('station_id'):
        r = E.change(g, boot=False)
        if r: out[sid] = r['change_pct']
    s = pd.Series(out).sort_values(); nxt = s.index[1] if s.index[0] == 'site_273' else s.index[0]
    L.append(f"| {yr or 'none'} | {len(s)} | {fm(s['site_273'])}% | {int(s.rank()['site_273'])} | {st.station_name[nxt]} {fm(s[nxt])}% |")
# raw station PM2.5 rank
out = {}
for sid, g in mp.groupby('station_id'):
    g = g.copy(); g['lr'] = g.l_st; out[sid] = E.change(g, boot=False)['change_pct']
s = pd.Series(out).sort_values()
L += ['', f"Station PM2.5 on its own: Ardhali Bazar {fm(s['site_273'])}%, rank {int(s.rank()['site_273'])} of {len(s)}; next {st.station_name[s.index[1]]} {fm(s.iloc[1])}%."]
open(f'{H}/gap_by_year.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
