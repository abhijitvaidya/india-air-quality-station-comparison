#!/usr/bin/env python3
"""PM2.5, PM10 and the ratio of PM2.5 to PM10 on the same station-months, 2023-2024 against 2018-2019 (Section 4.2,
Table A4). Months in which PM2.5 exceeds PM10 by more than 5% are left out.
Uses the estimator and data of 20_pm10_manual_checks.py.

Run:  python code/22_pm25_pm10_ratio.py
Output: results/pm25_pm10_ratio.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import numpy as np, pandas as pd
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '20_pm10_manual_checks.py')).read()
exec(src[:src.index("pc = lambda v")].replace("yv = g.lr.values - g.groupby(cell).lr.transform('mean').values", "yv = g[COL].values - g.groupby(cell)[COL].transform('mean').values"))
p = pd.read_csv(f'{DS}/continuous_station_month_pm25.csv', dtype={'ym': str}); p = p[(p.ndays >= 15) & (p['mean'] > 0)][['station_id', 'ym', 'mean']].rename(columns={'mean': 'pm25'})
b = c.merge(p, on=['station_id', 'ym']); b['l25'] = np.log(b.pm25); b['l10'] = np.log(b['mean']); b['frac'] = b.l25 - b.l10
b = b[b.pm25 <= b['mean'] * 1.05]          # PM2.5 cannot exceed PM10; drop the few months where it does by more than 5%
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']; upb = b.state == 'Uttar Pradesh'
G = {'Uttar Pradesh': b[upb], 'Rest of India': b[~upb], 'Plain outside Uttar Pradesh': b[b.state.isin(IGP)]}
pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
yrs = list(range(2018, 2025)); base = [2018, 2019]; i = [yrs.index(2023), yrs.index(2024)]
L = ['# PM2.5 and PM10 on the same station-months, 2023–2024 against 2018–2019', '', '| Group | Stations | Median PM2.5/PM10 in 2018–2019 | Median PM2.5/PM10 in 2023–2024 | PM2.5 | PM10 | PM2.5/PM10 ratio (interval) |', '|---|---|---|---|---|---|---|']
for name, g in G.items():
    out = []
    for COL in ['l25', 'l10', 'frac']:
        e, dr, _, ns, nk = estimate(g, yrs, base); out.append((pc(e[i].mean()), dr[:, i].mean(1)))
    v = out[2][1]; m0 = np.exp(g[g.y.isin(base)].frac).median(); m1 = np.exp(g[g.y.isin([2023, 2024])].frac).median()
    L.append(f'| {name} | {ns} | {m0:.2f} | {m1:.2f} | {fm(out[0][0])}% | {fm(out[1][0])}% | {fm(out[2][0])}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}) |')
L += ['', 'By year, PM2.5/PM10 ratio relative to 2018–2019:', '', '| Year | Uttar Pradesh | Rest of India | Plain outside UP |', '|---|---|---|---|']
COL = 'frac'; E = {k: estimate(g, yrs, base)[0] for k, g in G.items()}
for j, y in enumerate(yrs): L.append(f'| {y} | ' + ' | '.join(fm(pc(E[k][j])) + '%' for k in G) + ' |')
open(f'{RESULTS}/pm25_pm10_ratio.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
