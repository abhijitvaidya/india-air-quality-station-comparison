"""Two-period change at one station: season-balanced mean log ratio, 2016-2019 against 2023-2024, with a
season-year block bootstrap. Used by the single-station scripts."""
import math
import numpy as np
from common import SEASONS
B = 2000
rng = np.random.default_rng(20261117)


def tquant(p, df):
    # inverse Student t by bisection on the CDF (regularised incomplete beta via continued fraction)
    def betacf(a, b, x):
        qab, qap, qam = a + b, a + 1, a - 1; c, d = 1.0, 1 - qab * x / qap; d = 1 / d; h = d
        for mm in range(1, 300):
            m2 = 2 * mm; aa = mm * (b - mm) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d; c = 1 + aa / c; d = 1 / d; h *= d * c
            aa = -(a + mm) * (qab + mm) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d; c = 1 + aa / c; d = 1 / d; de = d * c; h *= de
            if abs(de - 1) < 1e-12: break
        return h
    def ibeta(a, b, x):
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
        return bt * betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1 - bt * betacf(b, a, 1 - x) / b
    def cdf(t):
        x = df / (df + t * t); v = 0.5 * ibeta(df / 2, 0.5, x)
        return 1 - v if t > 0 else v
    lo, hi = -1e3, 1e3
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(mid) < p: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def season_year(g, winter):
    if winter == 'calendar':
        return g.y.values
    return np.where(g.mo.isin([11, 12]), g.y + 1, g.y)


def draws(g, common, sy, kind='block_corr'):
    cols = []
    for s in common:
        msk = (g.season == s).values; v = g.lr.values[msk]; syv = sy[msk]
        point = v.mean()
        ub = np.unique(syv); nb = len(ub)
        if nb == 1 or kind == 'iid':
            cols.append(v[rng.integers(0, len(v), size=(B, len(v)))].mean(1)); continue
        bs = np.array([v[syv == u].sum() for u in ub]); bc = np.array([(syv == u).sum() for u in ub])
        idx = rng.integers(0, nb, size=(B, nb))
        dr = bs[idx].sum(1) / bc[idx].sum(1)
        if kind == 'block_corr':
            dr = point + (dr - point) * math.sqrt(nb / (nb - 1))
        cols.append(dr)
    return np.column_stack(cols).mean(1)


def change(g, pre=(2016, 2019), post=(2023, 2024), winter='following', drop_partial=False, min_n=2, need4=False,
           common_months=False, boot=True, kind='block_corr'):
    sy = season_year(g, winter)
    a = (g.y >= pre[0]) & (g.y <= pre[1]); b = (g.y >= post[0]) & (g.y <= post[1])
    if drop_partial and winter == 'following':
        # drop winter blocks that are incomplete by construction at the period edges
        part = ((g.season == 'winter') & (((sy == pre[0]) & a) | ((sy == pre[1] + 1) & a) | ((sy == post[0]) & b) | ((sy == post[1] + 1) & b))).values
        a = a & ~part; b = b & ~part
    if common_months:
        ma, mb = set(g[a].mo), set(g[b].mo); cm = ma & mb
        a = a & g.mo.isin(cm); b = b & g.mo.isin(cm)
    ga, gb = g[a], g[b]
    if len(ga) < 6 or len(gb) < 6:
        return None
    common = [s for s in SEASONS if (ga.season == s).sum() >= min_n and (gb.season == s).sum() >= min_n]
    if len(common) < (4 if need4 else 2):
        return None
    pa = np.mean([ga[ga.season == s].lr.mean() for s in common]); pb = np.mean([gb[gb.season == s].lr.mean() for s in common])
    out = {'months_ref': len(ga), 'months_late': len(gb), 'seasons': len(common), 'ratio_pre': math.exp(pa), 'ratio_post': math.exp(pb),
           'change_pct': 100 * (math.exp(pb - pa) - 1)}
    sya, syb = sy[a.values], sy[b.values]
    out['blocks_pre'] = ';'.join(f"{s[:3]}{len(np.unique(sya[(ga.season == s).values]))}" for s in common)
    out['blocks_post'] = ';'.join(f"{s[:3]}{len(np.unique(syb[(gb.season == s).values]))}" for s in common)
    out['min_blocks_post'] = min(len(np.unique(syb[(gb.season == s).values])) for s in common)
    out['fallback_strata'] = sum(len(np.unique(sya[(ga.season == s).values])) == 1 for s in common) + sum(len(np.unique(syb[(gb.season == s).values])) == 1 for s in common)
    if boot:
        d = draws(gb, common, syb, kind) - draws(ga, common, sya, kind)
        out['lo'] = 100 * (math.exp(np.quantile(d, 0.025)) - 1); out['hi'] = 100 * (math.exp(np.quantile(d, 0.975)) - 1)
    # year-level interval: season-balanced annual log ratios, Welch t
    def annual(gg, years):
        vals = []
        for yy in years:
            x = gg[gg.y == yy]
            ss = [x[x.season == s].lr.mean() for s in common if (x.season == s).sum() >= 1]
            if len(ss) == len(common): vals.append(np.mean(ss))
        return np.array(vals)
    ya = annual(ga, range(pre[0], pre[1] + 1)); yb = annual(gb, range(post[0], post[1] + 1))
    out['years_pre_complete'] = len(ya); out['years_post_complete'] = len(yb)
    if len(ya) >= 2 and len(yb) >= 2:
        va, vb = ya.var(ddof=1) / len(ya), yb.var(ddof=1) / len(yb); se = math.sqrt(va + vb)
        df = (va + vb) ** 2 / (va ** 2 / (len(ya) - 1) + vb ** 2 / (len(yb) - 1)) if va + vb > 0 else 1
        q = tquant(0.975, max(df, 1.0)); dm = yb.mean() - ya.mean()
        out['year_t_lo'] = 100 * (math.exp(dm - q * se) - 1); out['year_t_hi'] = 100 * (math.exp(dm + q * se) - 1); out['year_t_df'] = df
    return out


