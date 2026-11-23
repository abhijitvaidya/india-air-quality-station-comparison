#!/usr/bin/env python3
"""Two-period change at single stations, 2016-2019 to 2023-2024, against satellite PM2.5, CAMS PM2.5 and aerosol optical
depth (Section 3.3, Table F1, Figure 2).

Season-balanced mean log ratio (two_period.change) with a season-year block bootstrap. Intervals are widened on the
log scale by 1.62, the factor from the coverage simulation (26_interval_coverage.py). A station is eligible with at least
six valid months in each period and at least two seasons with two or more months in both periods.

Run:  python code/13_single_station.py
Outputs: results/single_station.csv, results/single_station.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys
import numpy as np, pandas as pd, netCDF4 as nc
H = str(RESULTS)
from common import paired
import two_period as E
W = 1.62
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
m['washu'] = m.sat_pm25
m['cams'] = pd.Series([cams.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float)
m['ceres'] = pd.Series([aod.get(k) for k in zip(m.station_id, m.ym)], index=m.index, dtype=float)
rows = []
for sid, g0 in m.groupby('station_id'):
    row = {'station_id': sid, 'station_name': g0.station_name.iloc[0], 'city': g0.city.iloc[0], 'state': g0.state.iloc[0], 'coord_source': st.coord_source[sid]}
    for ref in ['washu', 'cams', 'ceres']:
        g = g0.dropna(subset=[ref]).copy(); g = g[g[ref] > 0]; g['lr'] = np.log(g.cpcb / g[ref])
        r = E.change(g) if len(g) else None
        if r is None: continue
        c = np.log1p(r['change_pct'] / 100); lo = c + W * (np.log1p(r['lo'] / 100) - c); hi = c + W * (np.log1p(r['hi'] / 100) - c)
        row.update({f'{ref}_chg': r['change_pct'], f'{ref}_lo': 100 * np.expm1(lo), f'{ref}_hi': 100 * np.expm1(hi), f'{ref}_seasons': r['seasons']})
    if 'washu_chg' in row: rows.append(row)
x = pd.DataFrame(rows)
for ref in ['washu', 'cams', 'ceres']: x[f'{ref}_rank'] = x[f'{ref}_chg'].rank(method='min')
x.sort_values('washu_chg').round(2).to_csv(f'{H}/single_station.csv', index=False)
fm = lambda v: f'{v:+.0f}%'.replace('-', '−')
L = ['# Two-period change at single stations, 2016–2019 to 2023–2024', '',
     f'{len(x)} stations qualify against satellite PM2.5, {x.cams_chg.notna().sum()} against CAMS, {x.ceres_chg.notna().sum()} against optical depth. '
     f'Uttar Pradesh: {(x.state == "Uttar Pradesh").sum()}. Intervals widened by {W}.', '',
     '## Ten largest declines against satellite PM2.5', '', '| Station | State | Satellite PM2.5 (rank) | CAMS (rank) | Optical depth (rank) |', '|---|---|---|---|---|']
cell = lambda r, ref: 'n/a' if pd.isna(r[f'{ref}_chg']) else f"{fm(r[f'{ref}_chg'])} ({fm(r[f'{ref}_lo'])} to {fm(r[f'{ref}_hi'])}), {int(r[f'{ref}_rank'])}"
for _, r in x.sort_values('washu_chg').head(10).iterrows():
    L.append(f"| {r.station_name} | {r.state} | {cell(r, 'washu')} | {cell(r, 'cams')} | {cell(r, 'ceres')} |")
up = x.state == 'Uttar Pradesh'
L += ['', '## Medians', '', '| Group | Stations | Satellite PM2.5 | CAMS | Optical depth |', '|---|---|---|---|---|']
for nme, q in [('Uttar Pradesh', x[up]), ('Elsewhere', x[~up])]:
    L.append(f"| {nme} | {len(q)} | {fm(q.washu_chg.median())} | {fm(q.cams_chg.median())} | {fm(q.ceres_chg.median())} |")
L += ['', f'Intervals wholly below zero: satellite PM2.5 {int((x.washu_hi < 0).sum())}, CAMS {int((x.cams_hi < 0).sum())}, optical depth {int((x.ceres_hi < 0).sum())}; '
          f'wholly above: satellite PM2.5 {int((x.washu_lo > 0).sum())}, CAMS {int((x.cams_lo > 0).sum())}, optical depth {int((x.ceres_lo > 0).sum())}.']
open(f'{H}/single_station.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
