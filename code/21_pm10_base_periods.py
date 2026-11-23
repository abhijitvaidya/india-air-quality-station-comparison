#!/usr/bin/env python3
"""Continuous PM10 against CAMS PM10 with 2016-2019 and with 2018-2019 as the base (Section 4.2, Table A3).
Uses the estimator and data of 20_pm10_manual_checks.py.

Run:  python code/21_pm10_base_periods.py
Output: results/pm10_base_periods.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, re, numpy as np, pandas as pd
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '20_pm10_manual_checks.py')).read()
exec(src[:src.index("pc = lambda v")].replace("yv = g.lr.values - g.groupby(cell).lr.transform('mean').values", "yv = g[COL].values - g.groupby(cell)[COL].transform('mean').values"))
c['l_st'] = np.log(c['mean']); c['l_ref'] = c.l_st - c.lr
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']
G = {'All India': c, 'Uttar Pradesh': c[up], 'Rest of India': c[~up], 'Plain outside Uttar Pradesh': c[c.state.isin(IGP)]}
pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
Y = list(range(2016, 2025)); L = ['# Continuous PM10 against CAMS PM10, mean of 2023 and 2024, by base period', '']
for base, yrs in [([2016, 2017, 2018, 2019], Y), ([2018, 2019], Y[2:])]:
    L += [f'## Base {base[0]}–{base[-1]}', '', '| Group | Stations | Stations in the first base year | Station PM10 | CAMS PM10 at the stations | Station against CAMS (interval) |', '|---|---|---|---|---|---|']
    D = {}; i = [yrs.index(2023), yrs.index(2024)]
    for name, g in G.items():
        out = []
        for COL in ['l_st', 'l_ref', 'lr']:
            e, dr, _, ns, nk = estimate(g, yrs, base); out.append((pc(e[i].mean()), dr[:, i].mean(1)))
        D[name] = out[2]; v = out[2][1]
        L.append(f'| {name} | {ns} | {g[g.y == base[0]].station_id.nunique()} | {fm(out[0][0])}% | {fm(out[1][0])}% | {fm(out[2][0])}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}) |')
    for other in ['Rest of India', 'Plain outside Uttar Pradesh']:
        dd = pc(D['Uttar Pradesh'][1]) - pc(D[other][1])
        L.append(f'\nUttar Pradesh minus {other}: {fm(D["Uttar Pradesh"][0] - D[other][0])} points ({fm(np.nanpercentile(dd, 2.5))} to {fm(np.nanpercentile(dd, 97.5))}).')
    L.append('')
open(f'{RESULTS}/pm10_base_periods.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
