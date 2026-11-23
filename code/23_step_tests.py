#!/usr/bin/env python3
"""Single step at an unknown date in each station's monthly log ratio to satellite PM2.5 (Section 3.4, Appendix A6).
Stations with at least 60 paired months. Supremum F statistic (Andrews, 1993); reference distribution from 999 series
built from the fit without a step plus resampled residuals of the step model, in circular blocks. Detections at a 5% false
discovery rate (Benjamini-Hochberg, BH; Benjamini-Yekutieli, BY).

Specifications (one change at a time from 'orig'):
  block_len   residual block length in retained observations (orig 6)
  trim        minimum segment length in observations (orig 6)
  scheme      'rows' blocks over retained observations (orig); 'calendar' blocks over the station's calendar span;
              'sync' calendar blocks with the same block dates at every station
  null        'season' (orig) or 'trend' (a linear trend in both models; step beyond the trend)
  drop        calendar months excluded (2020; March to June 2020)
  covid       a March to June 2020 term in both models
  resid       'alt' residuals from the step fit (orig) or 'null' from the fit without a step
  seas        four season terms (orig) or twelve month terms

Run:  python code/23_step_tests.py <specification> [v6|v5]   (see run_all.sh for the list)
Outputs: results/step_tests/break_<specification>_<version>.csv, and the bootstrap draws for orig and sync
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json, os, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import paired, bh, by, SEASONS, HERE

NB = 999
MIN_MONTHS = 60
LOCK = ['202003', '202004', '202005', '202006']
VARIANTS = {
    'orig': {},
    'L3': {'block_len': 3}, 'L12': {'block_len': 12}, 'L24': {'block_len': 24},
    'trim12': {'trim': 12}, 'trim24': {'trim': 24},
    'calendar': {'scheme': 'calendar'}, 'calendar_L12': {'scheme': 'calendar', 'block_len': 12},
    'sync': {'scheme': 'sync'},
    'trend': {'null': 'trend'}, 'trend_calendar_L12_trim12': {'null': 'trend', 'scheme': 'calendar', 'block_len': 12, 'trim': 12},
    'drop_lockdown': {'drop': LOCK}, 'drop_2020': {'drop': [f'2020{m:02d}' for m in range(1, 13)]},
    'covid_dummy': {'covid': True},
    'resid_null': {'resid': 'null'},
    'combined': {'scheme': 'calendar', 'block_len': 12, 'trim': 12, 'drop': LOCK},
    'month12': {'seas': 'month12'}, 'month12_calendar_L12': {'seas': 'month12', 'scheme': 'calendar', 'block_len': 12},
    'month12_calendar_L6': {'seas': 'month12', 'scheme': 'calendar'},
    'combined_month12': {'seas': 'month12', 'scheme': 'calendar', 'block_len': 12, 'trim': 12, 'drop': LOCK},
    'combined_trend': {'scheme': 'calendar', 'block_len': 12, 'trim': 12, 'drop': LOCK, 'null': 'trend'},
}


def run(name, satv='v6', save_draws=False):
    o = dict(block_len=6, trim=6, scheme='rows', null='season', drop=[], covid=False, resid='alt', seas='season4')
    o.update(VARIANTS[name])
    rng = np.random.default_rng(20261025)
    TRIM, L = o['trim'], o['block_len']
    m = paired(satv)
    elig = set(m.groupby('station_id').size().loc[lambda x: x >= MIN_MONTHS].index)  # eligibility before any drop
    if o['drop']:
        m = m[~m.ym.isin(o['drop'])]

    def design0(g):
        if o['seas'] == 'month12':
            S = np.column_stack([(g.mo == k).values.astype(float) for k in range(1, 13)])
        else:
            S = np.column_stack([(g.season == s).values.astype(float) for s in SEASONS])
        cols = [S]
        if o['null'] == 'trend':
            tt = (g.t.values - g.t.values.mean()) / 120.0
            cols.append(tt[:, None])
        if o['covid']:
            cv = g.ym.isin(LOCK).values.astype(float)
            if cv.sum() > 0:
                cols.append(cv[:, None])
        return np.column_stack(cols)

    def sup_break(y, X0):
        # generalised version of the original: projects out X0 (season dummies [+ trend, covid]) exactly
        n = len(y)
        P = np.linalg.pinv(X0)
        e0 = y - X0 @ (P @ y)
        rss0 = e0 @ e0
        ks = np.arange(TRIM, n - TRIM + 1)
        D = (np.arange(n)[None, :] >= ks[:, None]).astype(float)
        Dt = D - (D @ P.T) @ X0.T
        num = (Dt @ e0) ** 2; den = (Dt * Dt).sum(1)
        den[den < 1e-12] = np.nan
        rss = rss0 - num / den
        j = int(np.nanargmin(rss))
        return ks[j], ks, rss, rss0

    def fstat(rss0, rssmin, n, p):
        return (rss0 - rssmin) / (rssmin / (n - p))

    def block_resample(e, L):  # original
        n = len(e); nb = int(np.ceil(n / L))
        starts = rng.integers(0, n, size=nb)
        return np.concatenate([np.take(e, np.arange(s, s + L), mode='wrap') for s in starts])[:n]

    def cal_resample(R, L, starts=None):
        T = len(R); nb = int(np.ceil(T / L))
        if starts is None:
            starts = rng.integers(0, T, size=nb)
        return np.concatenate([np.take(R, np.arange(s, s + L), mode='wrap') for s in starts])[:T]

    groups = [(sid, g.sort_values('t')) for sid, g in m.groupby('station_id')]
    groups = [(sid, g) for sid, g in groups if sid in elig]
    prep = []
    for sid, g in groups:
        n = len(g); y = g.lr.values
        X0 = design0(g)
        p0 = np.linalg.matrix_rank(X0)
        k, ks, rss, rss0 = sup_break(y, X0)
        F = fstat(rss0, rss[ks == k][0], n, p0 + 1)
        X = np.column_stack([X0, (np.arange(n) >= k).astype(float)])
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        null_fit = X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]
        if o['resid'] == 'null':
            resid = y - null_fit
        prep.append(dict(sid=sid, g=g, n=n, y=y, X0=X0, p0=p0, k=k, F=F, beta=beta, resid=resid, null_fit=null_fit))

    T0 = 108
    Fb = {}
    if o['scheme'] == 'sync':
        starts_all = [rng.integers(0, T0, size=int(np.ceil(T0 / L))) for _ in range(NB)]
    for P_ in prep:
        g = P_['g']; F_blk = np.full(NB, np.nan)
        if o['scheme'] == 'rows':
            for b in range(NB):
                yb = P_['null_fit'] + block_resample(P_['resid'], L)
                _, _, rr, r0 = sup_break(yb, P_['X0'])
                F_blk[b] = fstat(r0, np.nanmin(rr), P_['n'], P_['p0'] + 1)
        else:
            # season (and trend/covid) mean for every calendar month of the grid
            if o['scheme'] == 'calendar':
                t_lo, t_hi = int(g.t.min()), int(g.t.max())
            else:
                t_lo, t_hi = 0, T0 - 1
            grid = np.arange(t_lo, t_hi + 1)
            R = np.full(len(grid), np.nan); R[g.t.values - t_lo] = P_['resid']
            # null-model mean at each grid month: season effect (+ trend) from H0 fit, covid term ignored outside
            mo = grid % 12 + 1
            seas = np.select([np.isin(mo, [11, 12, 1, 2]), np.isin(mo, [3, 4, 5]), np.isin(mo, [6, 7, 8, 9])], SEASONS[:3], 'post-monsoon')
            gg = pd.DataFrame({'t': grid, 'season': seas, 'mo': mo, 'ym': [f'{2016 + x // 12}{x % 12 + 1:02d}' for x in grid]})
            X0g = design0(gg)
            coef = np.linalg.lstsq(P_['X0'], P_['y'], rcond=None)[0]
            if X0g.shape[1] != len(coef):  # covid column absent on grid subset (cannot happen on full grid)
                X0g = X0g[:, :len(coef)]
            mu = X0g @ coef
            if o['drop']:
                dropmask = gg.ym.isin(o['drop']).values
            else:
                dropmask = np.zeros(len(grid), bool)
            for b in range(NB):
                st = starts_all[b] if o['scheme'] == 'sync' else None
                Rb = cal_resample(R, L, st)
                ok = ~np.isnan(Rb) & ~dropmask
                if ok.sum() < 2 * TRIM + 6:
                    continue
                yb = (mu + Rb)[ok]
                X0b = X0g[ok]
                _, _, rr, r0 = sup_break(yb, X0b)
                F_blk[b] = fstat(r0, np.nanmin(rr), ok.sum(), np.linalg.matrix_rank(X0b) + 1)
        Fb[P_['sid']] = F_blk
    rows = []
    for P_ in prep:
        g = P_['g']; fb = Fb[P_['sid']]; fb = fb[~np.isnan(fb)]
        first = g.iloc[0]
        rows.append({'station_id': P_['sid'], 'station_name': first.station_name, 'city': first.city, 'state': first.state,
                     'months': P_['n'], 'first': g.ym.iloc[0], 'last': g.ym.iloc[-1], 'break_month': g.ym.values[P_['k']],
                     'shift_pct': 100 * (np.exp(P_['beta'][-1]) - 1), 'supF': P_['F'],
                     'p': (1 + np.sum(fb >= P_['F'])) / (len(fb) + 1), 'nb_valid': len(fb)})
    r = pd.DataFrame(rows)
    r['bh_5pct'] = bh(r.p.values)
    r['by_5pct'] = by(r.p.values)
    if save_draws:
        np.save(f"{STEP}/break_draws_{name}_{satv}.npy", np.vstack([Fb[s] for s in r.station_id]))
    return r, o


if __name__ == '__main__':
    name = sys.argv[1]; satv = sys.argv[2] if len(sys.argv) > 2 else 'v6'
    t0 = time.time()
    r, o = run(name, satv, save_draws=(name in ('sync', 'orig')))
    r.to_csv(f"{STEP}/break_{name}_{satv}.csv", index=False)
    print(name, satv, f'{time.time() - t0:.0f}s', len(r), int(r.bh_5pct.sum()))
