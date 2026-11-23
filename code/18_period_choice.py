#!/usr/bin/env python3
"""Single-station change against satellite PM2.5 for other choices of base and end period (Table A2).

Run:  python code/18_period_choice.py
Outputs: results/period_choice.csv, results/period_choice.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, sys, numpy as np, pandas as pd
H = str(RESULTS)
from common import paired
import two_period as E
m = paired('v6')
W = [((2016, 2019), (2023, 2024)), ((2017, 2019), (2023, 2024)), ((2016, 2018), (2023, 2024)), ((2016, 2018), (2022, 2024)),
     ((2016, 2019), (2022, 2024)), ((2017, 2020), (2023, 2024)), ((2016, 2019), (2020, 2021)), ((2016, 2019), (2021, 2022)), ((2016, 2019), (2022, 2023))]
rows = []
for pre, post in W:
    res = {}
    for sid, g in m.groupby('station_id'):
        r = E.change(g, pre=pre, post=post, boot=False)
        if r: res[sid] = (r['change_pct'], g.state.iloc[0])
    s = pd.DataFrame(res, index=['chg', 'state']).T; s.chg = s.chg.astype(float); up = s.state == 'Uttar Pradesh'
    rows.append({'base': f'{pre[0]}–{pre[1]}', 'end': f'{post[0]}–{post[1]}', 'stations': len(s), 'up': int(up.sum()),
                 'median_up': s[up].chg.median(), 'median_rest': s[~up].chg.median(), 'ardhali': s.chg['site_273'], 'ardhali_rank': int(s.chg.rank()['site_273'])})
x = pd.DataFrame(rows); x.round(2).to_csv(f'{H}/period_choice.csv', index=False)
fm = lambda v: f'{v:+.0f}%'.replace('-', '−')
L = ['| Base period | End period | Stations (Uttar Pradesh) | Median change, Uttar Pradesh | Median change, elsewhere | Ardhali Bazar (rank) |', '|---|---|---|---|---|---|']
L += [f'| {r.base} | {r.end} | {r.stations} ({r.up}) | {fm(r.median_up)} | {fm(r.median_rest)} | {fm(r.ardhali)} ({r.ardhali_rank}) |' for r in x.itertuples()]
open(f'{H}/period_choice.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
