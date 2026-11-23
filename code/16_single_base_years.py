#!/usr/bin/env python3
"""Each single year from 2016 to 2019 as the base: station against reference in 2023-2024, PM2.5 against satellite PM2.5
and CAMS, and PM10 against CAMS (Table A3).

Run:  python code/16_single_base_years.py
Output: results/single_base_years.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import numpy as np, pandas as pd
pc = lambda v: 100 * np.expm1(v); fm = lambda v: f'{v:+.0f}'.replace('-', '−')
IGP = ['Punjab', 'Haryana', 'Chandigarh', 'Delhi', 'Bihar', 'West Bengal']
L = ['# Each single year as the base: station against reference in 2023–2024', '',
     '| Pollutant, reference | Base year | UP stations in base year | Rest stations in base year | Uttar Pradesh | Rest of India | Plain outside UP | UP minus rest (interval) | UP minus plain (interval) |', '|---|---|---|---|---|---|---|---|---|']
def row(label, base, gets, nu, nr):
    D = {k: gets(k) for k in ['Uttar Pradesh', 'Rest of India', 'Plain outside Uttar Pradesh']}
    def diff(o):
        dd = pc(D['Uttar Pradesh'][1]) - pc(D[o][1]); return f"{fm(pc(D['Uttar Pradesh'][0]) - pc(D[o][0]))} ({fm(np.nanpercentile(dd, 2.5))} to {fm(np.nanpercentile(dd, 97.5))})"
    L.append(f"| {label} | {base} | {nu} | {nr} | {fm(pc(D['Uttar Pradesh'][0]))}% | {fm(pc(D['Rest of India'][0]))}% | {fm(pc(D['Plain outside Uttar Pradesh'][0]))}% | {diff('Rest of India')} | {diff('Plain outside Uttar Pradesh')} |")
# ---- PM2.5
ns = {'__file__': __file__}; src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '14_gap_by_year.py')).read(); exec(src[:src.index("pc = lambda v")], ns)
m = ns['m']; upm = m.state == 'Uttar Pradesh'; MK = {'Uttar Pradesh': upm, 'Rest of India': ~upm, 'Plain outside Uttar Pradesh': m.state.isin(IGP)}
for ref, rn in [('washu', 'PM2.5, satellite PM2.5'), ('cams', 'PM2.5, CAMS')]:
    for b in [2016, 2017, 2018, 2019]:
        ns['YEARS'] = list(range(b, 2025)); ns['BASE'] = [b]; ns['bi'] = [0]; ns['keep'] = [i for i, y in enumerate(ns['YEARS']) if y != 2019]
        if b == 2019: ns['REF'] = 2019
        i = [ns['YEARS'].index(2023), ns['YEARS'].index(2024)]
        def gets(k):
            g = m[MK[k] & m.y.isin(ns['YEARS'])].dropna(subset=[f'l_{ref}']).copy(); g['lr_'] = g.l_st - g[f'l_{ref}']
            e, dr, n = ns['estimate'](g, 'lr_'); return e[i].mean(), dr[:, i].mean(1)
        row(rn, b, gets, m[upm & (m.y == b)].station_id.nunique(), m[~upm & (m.y == b)].station_id.nunique())
# ---- PM10
n2 = {'__file__': __file__}; src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '20_pm10_manual_checks.py')).read(); exec(src[:src.index("pc = lambda v")], n2)
c = n2['c']; upc = c.state == 'Uttar Pradesh'; CK = {'Uttar Pradesh': upc, 'Rest of India': ~upc, 'Plain outside Uttar Pradesh': c.state.isin(IGP)}
for b in [2016, 2017, 2018, 2019]:
    yrs = list(range(b, 2025)); i = [yrs.index(2023), yrs.index(2024)]
    def gets(k):
        e, dr, _, _, _ = n2['estimate'](c[CK[k]], yrs, [b]); return e[i].mean(), dr[:, i].mean(1)
    row('PM10, CAMS', b, gets, c[upc & (c.y == b)].station_id.nunique(), c[~upc & (c.y == b)].station_id.nunique())
open(f'{RESULTS}/single_base_years.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
