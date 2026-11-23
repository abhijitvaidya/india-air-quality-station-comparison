#!/usr/bin/env python3
"""Continuous PM2.5 against each reference with 2016-2019 and with 2018-2019 as the base (Table A3).
Uses the estimator and data of 14_gap_by_year.py.

Run:  python code/15_pm25_base_periods.py
Output: results/pm25_base_periods.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import numpy as np, pandas as pd
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '14_gap_by_year.py')).read()
exec(src[:src.index("pc = lambda v")])
pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']
masks = {'All India': pd.Series(True, index=m.index), 'Uttar Pradesh': up, 'Rest of India': ~up, 'Plain outside Uttar Pradesh': m.state.isin(IGP)}
RN = {'washu': 'Satellite PM2.5', 'cams': 'CAMS PM2.5', 'ceres': 'Optical depth'}
L = ['# Continuous PM2.5 against three references, mean of 2023 and 2024, by base period', '']
for base in [[2016, 2017, 2018, 2019], [2018, 2019]]:
    YEARS = list(range(base[0], 2025)); BASE = base; bi = [YEARS.index(y) for y in BASE]; keep = [i for i, y in enumerate(YEARS) if y != REF]
    i = [YEARS.index(2023), YEARS.index(2024)]
    L += [f'## Base {base[0]}–{base[-1]}', '']
    for ref in RN:
        L += [f'### {RN[ref]}', '', '| Group | Stations | Stations in the first base year | Station PM2.5 | Reference at the stations | Station against reference (interval) |', '|---|---|---|---|---|---|']
        D = {}
        for name, mk in masks.items():
            g = m[mk & m.y.isin(YEARS)].dropna(subset=[f'l_{ref}']).copy(); g['lr_'] = g.l_st - g[f'l_{ref}']
            a = estimate(g, 'l_st'); b = estimate(g, f'l_{ref}'); r = estimate(g, 'lr_'); v = r[1][:, i].mean(1); D[name] = (pc(r[0][i].mean()), v)
            L.append(f'| {name} | {r[2]} | {g[g.y == base[0]].station_id.nunique()} | {fm(pc(a[0][i].mean()))}% | {fm(pc(b[0][i].mean()))}% | {fm(D[name][0])}% ({fm(pc(np.nanpercentile(v, 2.5)))} to {fm(pc(np.nanpercentile(v, 97.5)))}) |')
        for other in ['Rest of India', 'Plain outside Uttar Pradesh']:
            dd = pc(D['Uttar Pradesh'][1]) - pc(D[other][1])
            L.append(f'\nUttar Pradesh minus {other}: {fm(D["Uttar Pradesh"][0] - D[other][0])} points ({fm(np.nanpercentile(dd, 2.5))} to {fm(np.nanpercentile(dd, 97.5))}).')
        L.append('')
open(f'{RESULTS}/pm25_base_periods.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
