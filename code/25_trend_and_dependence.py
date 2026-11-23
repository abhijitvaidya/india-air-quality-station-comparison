#!/usr/bin/env python3
"""Step against gradual trend by station (Bayesian information criterion; Table A7), correlation of step-model residuals by
distance (Table A8), and the count of detections when the same calendar blocks are resampled at every station
(Appendix A7).

Run after 23_step_tests.py:  python code/25_trend_and_dependence.py
Outputs: results/step_tests/trend_vs_step.csv, residual_correlation_by_distance.csv, step_test_dependence.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import paired, bh, SEASONS, HERE, DS
R = str(STEP)
rng = np.random.default_rng(20261115)
m = paired('v6')
orig = pd.read_csv(f'{R}/break_orig_v6.csv', dtype={'break_month': str})
ids = list(orig.station_id)
st = pd.read_csv(f'{DS}/stations.csv').drop_duplicates('station_id').set_index('station_id')
TRIM = 6

def bic(y, X):
    b = np.linalg.lstsq(X, y, rcond=None)[0]; e = y - X @ b; n = len(y)
    return n * np.log(e @ e / n) + X.shape[1] * np.log(n), b

rows = []; resid_grid = {}
for sid in ids:
    g = m[m.station_id == sid].sort_values('t'); y = g.lr.values; n = len(y)
    S = np.column_stack([(g.season == s).values.astype(float) for s in SEASONS])
    tt = (g.t.values - g.t.values.mean()) / 12.0   # years
    b_s, _ = bic(y, S)
    b_t, bt = bic(y, np.column_stack([S, tt]))
    best_step = (np.inf, None, None); best_ts = (np.inf, None, None)
    for k in range(TRIM, n - TRIM + 1):
        d = (np.arange(n) >= k).astype(float)
        v1, c1 = bic(y, np.column_stack([S, d])); v1 += np.log(n)   # extra penalty for the estimated break date
        v2, c2 = bic(y, np.column_stack([S, tt, d])); v2 += np.log(n)
        if v1 < best_step[0]: best_step = (v1, k, c1)
        if v2 < best_ts[0]: best_ts = (v2, k, c2)
    cands = {'season only': b_s, 'trend': b_t, 'step': best_step[0], 'trend+step': best_ts[0]}
    best = min(cands, key=cands.get)
    rows.append({'station_id': sid, 'state': st.loc[sid, 'state'], 'city': st.loc[sid, 'city'], 'months': n,
                 'trend_pct_per_year': 100 * (np.exp(bt[-1]) - 1),
                 'bic_season': b_s, 'bic_trend': b_t, 'bic_step': best_step[0], 'bic_trend_step': best_ts[0], 'best_bic': best,
                 'step_vs_trend_dbic': best_step[0] - b_t})
    # residuals of the main step model on the 2016-2024 grid (for spatial correlation)
    o = orig.set_index('station_id').loc[sid]
    k = int(np.where(g.ym.values == o.break_month)[0][0])
    X = np.column_stack([S, (np.arange(n) >= k).astype(float)])
    e = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    Rg = np.full(108, np.nan); Rg[g.t.values] = e; resid_grid[sid] = Rg
tv = pd.DataFrame(rows).merge(orig[['station_id', 'shift_pct', 'bh_5pct', 'break_month']], on='station_id')
tv.to_csv(f'{R}/trend_vs_step.csv', index=False)
up = tv.state.eq('Uttar Pradesh')
L = ['# Step versus gradual trend, and dependence between stations', '',
     '## Which model fits best (BIC; the break date counts as an extra parameter)', '',
     '| group | stations | season only | trend | step | trend+step | median trend %/yr |', '|---|---|---|---|---|---|---|']
for lab, msk in [('Uttar Pradesh', up), ('Outside UP', ~up), ('UP, downward detections in the main specification', up & tv.bh_5pct & (tv.shift_pct < 0)),
                 ('Outside UP, detections in the main specification', ~up & tv.bh_5pct)]:
    x = tv[msk]; vc = x.best_bic.value_counts()
    L.append(f"| {lab} | {len(x)} | {vc.get('season only', 0)} | {vc.get('trend', 0)} | {vc.get('step', 0)} | {vc.get('trend+step', 0)} | {x.trend_pct_per_year.median():+.1f} |")
a = tv[tv.station_id == 'site_273'].iloc[0]
L += ['', f"Ardhali Bazar: best model {a.best_bic}; trend-only slope {a.trend_pct_per_year:+.1f}%/yr; BIC(step) - BIC(trend) = {a.step_vs_trend_dbic:.1f} (negative favours a step).",
      f"UP stations with a negative trend-only slope: {int((up & (tv.trend_pct_per_year < 0)).sum())} of {int(up.sum())}; outside UP: {int((~up & (tv.trend_pct_per_year < 0)).sum())} of {int((~up).sum())}.", '']

# ---------------- spatial dependence: residual correlation by distance
ll = st.loc[ids, ['lat', 'lon']].values
def km(a, b):
    la1, lo1, la2, lo2 = map(np.radians, [a[0], a[1], b[0], b[1]])
    return 2 * 6371 * np.arcsin(np.sqrt(np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2))
pairs = []
for i in range(len(ids)):
    for j in range(i + 1, len(ids)):
        a_, b_ = resid_grid[ids[i]], resid_grid[ids[j]]; ok = ~np.isnan(a_) & ~np.isnan(b_)
        if ok.sum() >= 36:
            pairs.append({'i': ids[i], 'j': ids[j], 'km': km(ll[i], ll[j]), 'same_city': st.loc[ids[i], 'city'] == st.loc[ids[j], 'city'],
                          'r': np.corrcoef(a_[ok], b_[ok])[0, 1], 'overlap': int(ok.sum())})
pc = pd.DataFrame(pairs)
pc['band'] = pd.cut(pc.km, [-1, 25, 50, 150, 300, 600, 5000], labels=['0-25', '25-50', '50-150', '150-300', '300-600', '600+'])
pc.groupby('band', observed=True).r.agg(['size', 'median', 'mean']).to_csv(f'{R}/residual_correlation_by_distance.csv')
L += ['## Correlation of step-model residuals between station pairs (monthly, >= 36 shared months)', '',
      '| distance km | pairs | median r | mean r |', '|---|---|---|---|']
for b, x in pc.groupby('band', observed=True):
    L.append(f'| {b} | {len(x)} | {x.r.median():.2f} | {x.r.mean():.2f} |')
L.append(f"| same city (any distance) | {int(pc.same_city.sum())} | {pc[pc.same_city].r.median():.2f} | {pc[pc.same_city].r.mean():.2f} |")

# ---------------- synchronized null draws: dependence of detection counts under the global null
D = np.load(f'{R}/break_draws_sync_v6.npy')   # stations x draws, same block dates for every station
sync = pd.read_csv(f'{R}/break_sync_v6.csv')
assert list(sync.station_id) == ids
nst, nb = D.shape
P = np.full(D.shape, np.nan)
for i in range(nst):
    f = D[i]; ok = ~np.isnan(f); fv = f[ok]
    ranks = (fv[None, :] >= fv[:, None]).sum(1) - 1        # others at least as large
    P[i, ok] = (1 + ranks) / len(fv)
cnt05 = np.nansum(P <= 0.05, 0)
bhc = np.array([bh(np.nan_to_num(P[:, b], nan=1.0)).sum() for b in range(nb)])
p0 = 0.05; var_ind = nst * p0 * (1 - p0)
de = cnt05.var() / var_ind
L += ['', '## Global-null simulation with dependence kept (synchronized calendar blocks, 999 draws)', '',
      f'Each draw resamples the same calendar blocks for all {nst} stations, so cross-station correlation of the residuals is kept.',
      f'Stations with p <= 0.05 per draw: mean {cnt05.mean():.1f} (independence: {nst * p0:.1f}); SD {cnt05.std():.1f} (independence: {np.sqrt(var_ind):.1f}); 95th percentile {np.percentile(cnt05, 95):.0f}; max {cnt05.max()}.',
      f'Variance inflation of the count: {de:.1f}; implied effective number of independent stations for counting purposes: {nst / de:.0f} of {nst}.',
      f'BH 5% rejections per null draw: mean {bhc.mean():.2f}; share of draws with any rejection {np.mean(bhc > 0):.3f}; 99th percentile {np.percentile(bhc, 99):.0f}; max {bhc.max()}.',
      f'Observed (sync p-values): {int(sync.bh_5pct.sum())} BH rejections; {int((sync.p <= 0.05).sum())} with p <= 0.05.', '']

# ---------------- city-cluster view of the UP contrast
o = orig.copy(); o['up'] = o.state.eq('Uttar Pradesh')
o['score'] = np.where(o.bh_5pct & (o.shift_pct < 0), -1, np.where(o.bh_5pct & (o.shift_pct > 0), 1, 0))
cty = o.groupby(['state', 'city']).agg(n=('score', 'size'), down=('score', lambda s: int((s < 0).sum())), upc=('score', lambda s: int((s > 0).sum())),
                                       mean_step=('shift_pct', 'mean'), up=('up', 'first')).reset_index()
cu, co = cty[cty.up], cty[~cty.up]
L += ['## City clusters', '',
      f'UP: {len(cu)} cities with tested stations; {int((cu.down > 0).sum())} have at least one downward detection, {int((cu.upc > 0).sum())} an upward one.',
      f'Outside UP: {len(co)} cities; {int((co.down > 0).sum())} with a downward detection, {int((co.upc > 0).sum())} with an upward one.',
      'UP cities: ' + '; '.join(f"{r.city} {r.down}/{r.n}" for _, r in cu.sort_values('n', ascending=False).iterrows()), '']
# permutation over cities: state label shuffled among cities (city = unit), statistic = UP mean of city net detection share
cty['net'] = (cty.upc - cty.down) / cty.n
obs = cty[cty.up].net.mean() - cty[~cty.up].net.mean()
lab = cty.up.values.copy(); sims = []
for _ in range(20000):
    rng.shuffle(lab); sims.append(cty.net[lab].mean() - cty.net[~lab].mean())
sims = np.array(sims)
obs_step = cty[cty.up].mean_step.mean() - cty[~cty.up].mean_step.mean()
lab = cty.up.values.copy(); sims2 = []
for _ in range(20000):
    rng.shuffle(lab); sims2.append(cty.mean_step[lab].mean() - cty.mean_step[~lab].mean())
sims2 = np.array(sims2)
L += [f'City-level permutation (cities as units, {len(cty)} cities, UP label shuffled across cities, 20,000 draws; exploratory, the UP contrast was chosen after seeing the data):',
      f'- net detection share (up minus down, per city): UP mean {cty[cty.up].net.mean():+.2f}, others {cty[~cty.up].net.mean():+.2f}; one-sided p = {np.mean(sims <= obs):.4f}',
      f'- mean step estimate per city: UP {cty[cty.up].mean_step.mean():+.1f}%, others {cty[~cty.up].mean_step.mean():+.1f}%; one-sided p = {np.mean(sims2 <= obs_step):.4f}']
open(f'{R}/step_test_dependence.md', 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
