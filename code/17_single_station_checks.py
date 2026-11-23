#!/usr/bin/env python3
"""Checks on the single-station comparison (Appendix A2, A8; Table A1; Figure 3).

1. Two-period change against satellite PM2.5 under five other eligibility rules (Table A1).
2. Leave one year out: range of the median, sign changes, rank of Ardhali Bazar.
3. Season-year blocks behind the intervals.
4. Each station's ratio to satellite PM2.5 in the two periods, and the change in its own PM2.5 (Table F1).
5. Stations falling against optical depth and CAMS, and the number of optical-depth grid cells.
6. Ardhali Bazar's season-balanced ratio by calendar year (Figure 3).

Run after 13_single_station.py:  python code/17_single_station_checks.py
Outputs: results/two_period_variants.csv, results/single_station_ratios.csv, results/ardhali_yearly_ratio.csv,
results/single_station_checks.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys
import numpy as np, pandas as pd
H = str(RESULTS)
from common import paired, SEASONS
import two_period as E
m = paired('v6'); groups = {sid: g.sort_values('t') for sid, g in m.groupby('station_id')}
L = ['# Checks on the single-station comparison', '']
fm = lambda v, d=0: f'{v:+.{d}f}%'.replace('-', '−')

# 1. eligibility rules
VAR = {'orig': {}, 'drop_partial_winters': {'drop_partial': True}, 'winter_calendar_year': {'winter': 'calendar'},
       'all_four_seasons': {'need4': True}, 'min3_per_season': {'min_n': 3}, 'common_calendar_months': {'common_months': True}}
rows, base = [], None
for v, kw in VAR.items():
    d = pd.DataFrame([{'station_id': sid, 'state': g.state.iloc[0], **r} for sid, g in groups.items() if (r := E.change(g, boot=False, **kw))])
    if v == 'orig': base = d.set_index('station_id')
    a = d[d.station_id == 'site_273']
    rows.append({'variant': v, 'stations': len(d), 'up': int((d.state == 'Uttar Pradesh').sum()), 'median_change': d.change_pct.median(),
                 'ardhali_change': a.change_pct.iloc[0], 'ardhali_rank': int(d.change_pct.rank(method='min')[a.index].iloc[0])})
vt = pd.DataFrame(rows); vt.round(2).to_csv(f'{H}/two_period_variants.csv', index=False)
L += ['## 1. Two-period change under other eligibility rules', '', vt.round(1).to_string(index=False), '']

# 2. leave one year out
lo = []
for sid in base.index:
    g = groups[sid]
    for yy in [2016, 2017, 2018, 2019, 2023, 2024]:
        r = E.change(g[g.y != yy], boot=False); lo.append({'station_id': sid, 'drop': yy, 'chg': r['change_pct'] if r else np.nan})
lw = pd.DataFrame(lo).pivot(index='station_id', columns='drop', values='chg')
med = lw.median(); flips = int(((np.sign(lw).min(1) < 0) & (np.sign(lw).max(1) > 0)).sum())
rk = {yy: int(lw[yy].rank(method='min')['site_273']) for yy in lw.columns}
nxt = {yy: round(float(lw[yy].drop('site_273').min() - lw.loc['site_273', yy]), 1) for yy in lw.columns}
L += ['## 2. Leave one year out', '', f'Stations: {len(base)}. Median with all years {fm(base.change_pct.median(), 1)}; with one year left out, from {fm(med.max(), 1)} to {fm(med.min(), 1)}.',
      f'Stations whose change has a different sign in at least one of the six runs: {flips}.',
      f'Ardhali Bazar: {", ".join(f"{yy}: {fm(lw.loc[chr(115) + "ite_273", yy], 1)} (rank {rk[yy]}, next station {nxt[yy]} points behind)" for yy in lw.columns)}.', '']

# 3. blocks
b = base
a = b.loc['site_273']
L += ['## 3. Season-year blocks', '', f'Ardhali Bazar: 2016–2019 {a.blocks_pre}; 2023–2024 {a.blocks_post}.',
      f'Stations with at least two blocks in every compared season of 2023–2024: {int((b.min_blocks_post >= 2).sum())} of {len(b)}.',
      f'Seasons compared: {b.seasons.value_counts().sort_index().to_dict()}.', '']

# 4. ratios and the station's own change
cp = m.copy(); own = []
for sid in base.index:
    g = groups[sid].copy(); g['lr'] = np.log(g.cpcb); r1 = E.change(g, boot=False)
    g['lr'] = np.log(g.sat_pm25); r2 = E.change(g, boot=False)
    own.append({'station_id': sid, 'ratio_2016_19': base.ratio_pre[sid], 'ratio_2023_24': base.ratio_post[sid], 'station_chg': r1['change_pct'], 'satellite_chg': r2['change_pct']})
own = pd.DataFrame(own); own.round(4).to_csv(f'{H}/single_station_ratios.csv', index=False)
o = own.set_index('station_id'); ss = pd.read_csv(f'{H}/single_station.csv').set_index('station_id')
assert (o.index.sort_values() == ss.index.sort_values()).all()
L += ['## 4. Own change', '', f"Ardhali Bazar: ratio {o.ratio_2016_19['site_273']:.2f} to {o.ratio_2023_24['site_273']:.2f}; own PM2.5 {fm(o.station_chg['site_273'])} (rank {int(o.station_chg.rank(method='min')['site_273'])} of {len(o)}); "
      f"satellite estimate {fm(o.satellite_chg['site_273'])}.", '']

# 5. optical depth
up = ss.state == 'Uttar Pradesh'
st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id'); sig = set()
for sid in ss.index:
    f = aod_file(st.lat[sid], st.lon[sid])
    if os.path.exists(f):
        p = json.load(open(f))['properties']['parameter']['AOD_55']; sig.add(tuple(p[k] for k in sorted(p) if k[4:] != '13'))
L += ['## 5. Optical depth and CAMS', '', f'Against optical depth: {int((ss.ceres_chg[up] < 0).sum())} of {int(up.sum())} stations in Uttar Pradesh fall, {int((ss.ceres_chg[~up] < 0).sum())} of {int((~up).sum())} elsewhere.',
      f'Against CAMS: {int((ss.cams_chg[up] < 0).sum())} of {int(up.sum())}, {int((ss.cams_chg[~up] < 0).sum())} of {int((~up).sum())}.',
      f'Against the satellite estimate: {int((ss.washu_chg[up] < 0).sum())} of {int(up.sum())}, {int((ss.washu_chg[~up] < 0).sum())} of {int((~up).sum())}.',
      f'Distinct optical-depth series among the {len(ss)} stations: {len(sig)}.', '']

# 6. Ardhali Bazar, season-balanced ratio to satellite PM2.5 by calendar year (Figure 3)
g = groups['site_273']; ann = []
for yy in range(2016, 2025):
    x = g[g.y == yy]
    ann.append({'year': yy, 'ratio': float(np.exp(np.mean([x[x.season == s].lr.mean() for s in SEASONS]))), 'months': len(x)})
ann = pd.DataFrame(ann); ann.round(4).to_csv(f'{H}/ardhali_yearly_ratio.csv', index=False)
L += ['## 6. Ardhali Bazar by year', '', '; '.join(f'{int(r.year)}: {r.ratio:.2f} ({int(r.months)} months)' for r in ann.itertuples()), '']
open(f'{H}/single_station_checks.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
