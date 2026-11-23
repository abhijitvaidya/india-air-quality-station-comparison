#!/usr/bin/env python3
"""Coverage achieved by the widened interval (factor 1.62) when there is no true change, on fresh simulations
(Appendix A5, Table A5). Same data-generating model as 26_interval_coverage.py, on the same 159 stations; new seeds.

Run:  python code/27_achieved_coverage.py   (200 simulations a station; a few minutes)
Outputs: results/achieved_coverage.csv, results/achieved_coverage_station.csv, results/achieved_coverage.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = str(RESULTS)
import two_period as E          # change(), draws(); uses its own module-level rng
from common import paired

FACTOR = 1.62
NSIM = 200
E.rng = np.random.default_rng(20261118)       # fresh bootstrap draws
rng = np.random.default_rng(118118)           # fresh simulated series

m = paired('v6')
st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id')
m = m[m.station_id.map(st.coord_source) == 'database']   # the 159 stations of 26_interval_coverage.py
groups = {sid: g.sort_values('t') for sid, g in m.groupby('station_id')}
base = []
for sid, g in groups.items():
    r = E.change(g)
    if r:
        base.append({'station_id': sid, 'seasons': r['seasons'], 'fallback_strata': r['fallback_strata'], 'min_blocks_post': r['min_blocks_post']})
base = pd.DataFrame(base)

# parameters, as in 26_interval_coverage.py
pre = m[(m.y <= 2019) & m.station_id.isin(base.station_id)].copy()
pre['e'] = pre.lr - pre.groupby(['station_id', 'season']).lr.transform('mean')
ym = pre.groupby(['station_id', 'y']).e.agg(['mean', 'size']).reset_index()
pre = pre.merge(ym.rename(columns={'mean': 'ybar', 'size': 'ny'}), on=['station_id', 'y'])
pre['w'] = pre.e - pre.ybar
par = []
for sid, g in pre.groupby('station_id'):
    g = g.sort_values('t'); w = g.w.values; t = g.t.values
    adj = np.diff(t) == 1
    phi = np.corrcoef(w[:-1][adj], w[1:][adj])[0, 1] if adj.sum() > 10 else np.nan
    sw2 = w.var(ddof=1)
    yy = ym[ym.station_id == sid]
    sb2 = max(yy['mean'].var(ddof=1) - sw2 / yy['size'].mean(), 0) if len(yy) >= 3 else np.nan
    par.append({'station_id': sid, 'phi': phi, 'sw2': sw2, 'sb2': sb2})
par = pd.DataFrame(par)
PHI = float(np.nanmedian(par.phi)); SHARE = float(np.nanmedian(par.sb2 / (par.sb2 + par.sw2)))
pidx = par.set_index('station_id')

rows = []
for scen, share in [('estimated year share', SHARE), ('no year effect', 0.0), ('double year share', min(2 * SHARE, 0.9))]:
    for sid in base.station_id:
        g = groups[sid]; g = g[(g.y <= 2019) | (g.y >= 2023)].copy()
        tot = float(pidx.loc[sid, ['sw2', 'sb2']].fillna(0).sum()) if sid in pidx.index else float(np.nanmedian(par.sw2))
        sb2, sw2 = share * tot, (1 - share) * tot
        smean = g.groupby('season').lr.transform('mean').values
        tt = g.t.values; years = g.y.values
        raw = wid = 0
        for _ in range(NSIM):
            u = dict(zip(range(2016, 2025), rng.normal(0, math.sqrt(sb2), 9)))
            ar = np.zeros(108); z = rng.normal(0, math.sqrt(sw2 * (1 - PHI ** 2)), 108); ar[0] = rng.normal(0, math.sqrt(sw2))
            for i in range(1, 108):
                ar[i] = PHI * ar[i - 1] + z[i]
            gs = g.copy(); gs['lr'] = smean + np.array([u[yv] for yv in years]) + ar[tt]
            r = E.change(gs)
            c0 = math.log1p(r['change_pct'] / 100); lo = math.log1p(r['lo'] / 100); hi = math.log1p(r['hi'] / 100)
            raw += lo <= 0 <= hi
            wid += (c0 - FACTOR * (c0 - lo)) <= 0 <= (c0 + FACTOR * (hi - c0))
        rows.append({'scenario': scen, 'year_share': share, 'station_id': sid, 'raw': raw / NSIM, 'widened': wid / NSIM})
    print(scen, 'done', flush=True)
st = pd.DataFrame(rows).merge(base, on='station_id')
st.to_csv(f'{HERE}/achieved_coverage_station.csv', index=False)
summ = st.groupby('scenario', sort=False).agg(year_share=('year_share', 'first'), stations=('station_id', 'size'), raw=('raw', 'mean'), widened=('widened', 'mean')).reset_index()
summ.to_csv(f'{HERE}/achieved_coverage.csv', index=False)
L = ['# Achieved coverage of the widened interval', '',
     f'No true change; nominal 95%; {NSIM} fresh simulations for each of {len(base)} stations and each assumption. Widening factor {FACTOR} throughout.', '',
     '| Assumed year-to-year share of variance | Interval as run | Widened by 1.62 |', '|---|---|---|']
for r in summ.itertuples():
    L.append(f'| {r.scenario} ({r.year_share:.2f}) | {100 * r.raw:.1f}% | {100 * r.widened:.1f}% |')
e = st[st.scenario == 'estimated year share']
L += ['', '## Estimated year share, by how sparse the station is', '', '| Group | Stations | Interval as run | Widened by 1.62 |', '|---|---|---|---|']
for name, msk in [('Four seasons compared', e.seasons == 4), ('Three seasons', e.seasons == 3), ('Two seasons', e.seasons == 2),
                  ('At least one season with a single block', e.fallback_strata > 0), ('Every season with two or more blocks', e.fallback_strata == 0)]:
    L.append(f'| {name} | {int(msk.sum())} | {100 * e.raw[msk].mean():.1f}% | {100 * e.widened[msk].mean():.1f}% |')
L += ['', f"Lowest station coverage after widening (estimated share): {100 * e.widened.min():.0f}%; stations below 90%: {int((e.widened < 0.90).sum())} of {len(e)}.",
      f"Ardhali Bazar (site_273): as run {100 * e.raw[e.station_id == 'site_273'].iloc[0]:.0f}%, widened {100 * e.widened[e.station_id == 'site_273'].iloc[0]:.0f}%."]
open(f'{HERE}/achieved_coverage.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
