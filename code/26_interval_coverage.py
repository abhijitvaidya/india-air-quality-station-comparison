#!/usr/bin/env python3
"""Coverage of the two-period interval when there is no true change, and the widening factor (Appendix A5, Table A5).

Each eligible station keeps its own months and season means; the simulated log ratio adds a random effect for each year
and AR(1) monthly noise, with parameters estimated from 2016-2019 residuals. For each simulated series the nominal 95%
season-year block interval (two_period.change) is checked for zero. Three assumptions on the year-to-year share of
variance: as estimated, none, and twice the estimate. The widening factor is the one that would restore 95% under a
normal approximation: 1.96 / z, where z is the normal quantile of (1 + coverage) / 2.

The simulation uses the 159 eligible stations whose coordinates come from the station database; the factor is then
applied to all 166 eligible stations. A Welch t interval on yearly values is reported for comparison.

Run one assumption at a time, then without arguments to combine:
  python code/26_interval_coverage.py 0     (estimated share; 0 = estimated, 1 = none, 2 = twice)
  python code/26_interval_coverage.py 1
  python code/26_interval_coverage.py 2
  python code/26_interval_coverage.py
Outputs: results/interval_coverage_<assumption>.csv, results/interval_coverage.csv, results/interval_coverage.md
"""
import math, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
from common import paired
import two_period as E

NSIM = 100                       # simulations a station (27_achieved_coverage.py uses 200 fresh ones)
m = paired('v6')
st = pd.read_csv(f'{DERIVED}/stations.csv').set_index('station_id')
m = m[m.station_id.map(st.coord_source) == 'database']
groups = {sid: g.sort_values('t') for sid, g in m.groupby('station_id')}
base = pd.DataFrame([{'station_id': sid, **r} for sid, g in groups.items() if (r := E.change(g, boot=False))])

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
SCEN = [('estimated year share', SHARE), ('no year effect', 0.0), ('double year share', min(2 * SHARE, 0.9))]
fname = lambda scen: f"{RESULTS}/interval_coverage_{scen.replace(' ', '_')}.csv"

if len(sys.argv) > 1:
    k = int(sys.argv[1]); scen, share = SCEN[k]
    rng = np.random.default_rng(20261116 + 7 * k); E.rng = rng      # one stream for the series and the bootstrap
    hit = {'block_corr': [], 'block': [], 'iid': [], 'year_t': []}
    for sid in base.station_id:
        g = groups[sid]; g = g[(g.y <= 2019) | (g.y >= 2023)].copy()
        tot = float(par.set_index('station_id').loc[sid, ['sw2', 'sb2']].fillna(0).sum()) if sid in set(par.station_id) else float(np.nanmedian(par.sw2))
        sb2, sw2 = share * tot, (1 - share) * tot
        smean = g.groupby('season').lr.transform('mean').values
        tt = g.t.values; years = g.y.values
        for _ in range(NSIM):
            u = dict(zip(range(2016, 2025), rng.normal(0, math.sqrt(sb2), 9)))
            ar = np.zeros(108); z = rng.normal(0, math.sqrt(sw2 * (1 - PHI ** 2)), 108); ar[0] = rng.normal(0, math.sqrt(sw2))
            for i in range(1, 108): ar[i] = PHI * ar[i - 1] + z[i]
            gs = g.copy(); gs['lr'] = smean + np.array([u[yv] for yv in years]) + ar[tt]
            r_c = E.change(gs, kind='block_corr'); r_b = E.change(gs, kind='block', boot=True); r_i = E.change(gs, kind='iid')
            hit['block_corr'].append(r_c['lo'] <= 0 <= r_c['hi']); hit['block'].append(r_b['lo'] <= 0 <= r_b['hi'])
            hit['iid'].append(r_i['lo'] <= 0 <= r_i['hi'])
            if 'year_t_lo' in r_c: hit['year_t'].append(r_c['year_t_lo'] <= 0 <= r_c['year_t_hi'])
    pd.DataFrame([{'scenario': scen, 'year_share': share, 'phi': PHI, 'sims': NSIM * len(base), **{f'coverage_{k}': np.mean(v) for k, v in hit.items()}}]).to_csv(fname(scen), index=False)
    print(scen, 'done'); sys.exit(0)


def zinv(p):
    lo, hi = 0.0, 10.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if 0.5 * (1 + math.erf(mid / math.sqrt(2))) < p: lo = mid
        else: hi = mid
    return mid


cov = pd.concat([pd.read_csv(fname(s)) for s, _ in SCEN if os.path.exists(fname(s))])
cov['widening_factor'] = [1.96 / zinv((1 + c) / 2) for c in cov.coverage_block_corr]
cov.to_csv(f'{RESULTS}/interval_coverage.csv', index=False)
yt = base.dropna(subset=['year_t_lo']); a = base.set_index('station_id').loc['site_273']
L = ['# Coverage of the two-period interval with no true change', '',
     f'{len(base)} stations, {NSIM} simulations each. Year-to-year share of variance as estimated: {SHARE:.3f}; AR(1) coefficient {PHI:.3f}.', '',
     '| Assumption | Coverage of the nominal 95% interval | Widening factor | Welch t interval on yearly values |', '|---|---|---|---|']
L += [f'| {r.scenario} | {100 * r.coverage_block_corr:.1f}% | {r.widening_factor:.2f} | {100 * r.coverage_year_t:.1f}% |' for r in cov.itertuples()]
L += ['', f'Welch t interval exists for {len(yt)} of {len(base)} stations; median width {float((yt.year_t_hi - yt.year_t_lo).median()):.0f} points; '
          f'Ardhali Bazar {a.year_t_lo:.0f}% to {a.year_t_hi:.0f}%.']
open(f'{RESULTS}/interval_coverage.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
