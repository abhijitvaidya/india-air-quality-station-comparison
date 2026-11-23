#!/usr/bin/env python3
"""Weather adjustment of daily PM2.5 change from 2019 (Section 3.6, Appendix B, Table B1, Figure 4a).

log(daily mean PM2.5) is regressed on log(rain + 1), log(wind speed), temperature and relative humidity, with month,
weekday and year terms, over 2016-2024; the adjusted change compares the year terms with 2019. Intervals: 999 resamples
of whole calendar months within each year. Varanasi is Ardhali Bazar; the other cities are panels of the stations that
reported in both 2019 and 2024. A second fit on station-days with a term for each station removes the effect of a
changing panel.

Run:  python code/29_weather_adjustment.py   (needs data/raw/aq.duckdb and data/raw/weather; see data/README.md)
Outputs: results/weather_adjustment_results.csv, weather_adjustment_availability.csv, weather_adjustment_blh_gap.csv,
weather_adjustment_reconcile.csv, weather_adjustment.json
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json
import os

import duckdb
import numpy as np
import pandas as pd

OUT = str(RESULTS)
os.makedirs(OUT, exist_ok=True)
SEED = 20261115
B = 999
CITIES = ['Varanasi', 'Delhi', 'Patna', 'Lucknow', 'Kanpur']
rng = np.random.default_rng(SEED)

c = duckdb.connect(str(HOURLY_DB), read_only=True)
w_all = pd.read_csv(WEATHER_CSV, parse_dates=['d'])
PANEL_SQL = """with lng as (select station_id from daily_all where p='PM2.5' and source='cpcb_caaqm'
                 group by 1 having min(d) <= '2019-12-31' and max(d) >= '2024-06-30')"""
# station-day panel (original definition of the fixed panel)
sd = c.sql(PANEL_SQL + """ select city, station_id, d, v, n from daily_all join lng using (station_id)
  where p='PM2.5' and source='cpcb_caaqm' and d >= '2016-01-01'""").df()
sd['d'] = pd.to_datetime(sd.d)
city_day = sd.groupby(['city', 'd']).agg(pm=('v', 'median'), n_st=('v', 'size')).reset_index()

# ------------------------------------------------------------------ 1. BLH gap
gap = w_all[w_all.blh_m.isna()].groupby('city').d.agg(['min', 'max', 'count']).reset_index()
gap.to_csv(os.path.join(OUT, 'weather_adjustment_blh_gap.csv'), index=False)


def design(x, use_blh=True, weather=True, station_fe=False):
    base = 2019
    cols, names = [np.ones(len(x))], ['const']
    if weather:
        cols += [np.log1p(x.rain_mm.values), np.log(x.wind_kmh.values), x.temp_c.values, x.rh_pct.values]
        names += ['rain', 'wind', 'temp', 'rh']
        if use_blh:
            cols.append(np.log(x.blh_m.values)); names.append('blh')
    yrs = sorted(x.d.dt.year.unique())
    for y in yrs:
        if y != base:
            cols.append((x.d.dt.year == y).values.astype(float)); names.append(f'y{y}')
    for m in range(2, 13):
        cols.append((x.d.dt.month == m).values.astype(float)); names.append(f'm{m}')
    for dw in range(1, 7):
        cols.append((x.d.dt.dayofweek == dw).values.astype(float)); names.append(f'w{dw}')
    if station_fe:
        st = sorted(x.station_id.unique())
        for s in st[1:]:
            cols.append((x.station_id == s).values.astype(float)); names.append(f's_{s}')
    return np.column_stack(cols), names


def year_eff(X, names, yv):
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    return {int(n[1:]): float(b) for n, b in zip(names, beta) if n.startswith('y')}, beta


def boot(x, yv, use_blh, weather, station_fe, targets):
    """Resample calendar months within each year; refit; return draws of year effects for targets."""
    ym = (x.d.dt.year * 100 + x.d.dt.month).values
    yr = x.d.dt.year.values
    blocks = {}
    for y in np.unique(yr):
        mths = np.unique(ym[yr == y])
        blocks[y] = [np.where(ym == m)[0] for m in mths]
    draws = []
    for _ in range(B):
        idx = np.concatenate([np.concatenate([bl[j] for j in rng.integers(0, len(bl), len(bl))]) for bl in blocks.values()])
        xb = x.iloc[idx]
        Xb, nb = design(xb, use_blh, weather, station_fe)
        e, _ = year_eff(Xb, nb, yv[idx])
        draws.append([e.get(t, np.nan) for t in targets])
    return np.array(draws)


def pct(b):
    return 100 * (np.exp(b) - 1)


SPECS = {
    'A original (2016-2026, BLH, 2024 = Jan 1 + Jul-Dec)': dict(end='2026-12-31', blh=True, fe=False),
    'B main period 2016-2024, BLH (2024 = Jan 1 + Jul-Dec)': dict(end='2024-12-31', blh=True, fe=False),
    'C main period 2016-2024, no BLH (full 2024)': dict(end='2024-12-31', blh=False, fe=False),
    'D 2019-2024 only, no BLH': dict(start='2019-01-01', end='2024-12-31', blh=False, fe=False),
    'E station fixed effects, 2016-2024, no BLH': dict(end='2024-12-31', blh=False, fe=True),
}
TARGETS = [2023, 2024]
rows, recon, avail = [], [], []
for city in CITIES:
    w = w_all[w_all.city == city]
    for spec, s in SPECS.items():
        if s['fe']:
            x = sd[sd.city == city].merge(w, on=['city', 'd'])
        else:
            x = city_day[city_day.city == city].merge(w, on=['city', 'd'])
        x = x[(x.d >= s.get('start', '2016-01-01')) & (x.d <= s['end'])]
        need = ['rain_mm', 'wind_kmh', 'temp_c', 'rh_pct'] + (['blh_m'] if s['blh'] else [])
        x = x.dropna(subset=need)
        val = 'v' if s['fe'] else 'pm'
        x = x[(x[val] > 0) & (x.wind_kmh > 0) & ((x.blh_m > 0) if s['blh'] else True)].reset_index(drop=True)
        yv = np.log(x[val].values)
        out = {'city': city, 'spec': spec, 'obs': len(x), 'stations': x.station_id.nunique() if s['fe'] else None}
        for lab, weather in [('month_only', False), ('weather', True)]:
            X, nm = design(x, s['blh'], weather, s['fe'])
            e, beta = year_eff(X, nm, yv)
            dr = boot(x, yv, s['blh'], weather, s['fe'], TARGETS)
            for j, t in enumerate(TARGETS):
                out[f'{lab}_{t}'] = round(pct(e[t]), 1)
                lo, hi = np.nanpercentile(dr[:, j], [2.5, 97.5])
                out[f'{lab}_{t}_lo'], out[f'{lab}_{t}_hi'] = round(pct(lo), 1), round(pct(hi), 1)
            if weather:
                pred = X @ beta
                out['r2_with_year'] = round(1 - ((yv - pred) ** 2).sum() / ((yv - yv.mean()) ** 2).sum(), 2)
        # raw changes on the same observations (geometric and arithmetic)
        g = x.assign(y=x.d.dt.year, lv=np.log(x[val]))
        if not s['fe']:
            for t in TARGETS:
                if t in g.y.values and 2019 in g.y.values:
                    out[f'raw_geo_{t}'] = round(pct(g[g.y == t].lv.mean() - g[g.y == 2019].lv.mean()), 1)
                    out[f'raw_arith_{t}'] = round(100 * (g[g.y == t][val].mean() / g[g.y == 2019][val].mean() - 1), 1)
                    out[f'days_{t}'] = int((g.y == t).sum())
            out['days_2019'] = int((g.y == 2019).sum())
        rows.append(out)

    # reconciliation table, this city: annual all valid days vs model sample
    cd = city_day[city_day.city == city].copy()
    cd['y'] = cd.d.dt.year
    full = cd.merge(w, on=['city', 'd'], how='left')
    for y in [2019, 2023, 2024]:
        a = cd[cd.y == y]
        b = full[(full.y == y)].dropna(subset=['blh_m'])
        recon.append({'city': city, 'year': y, 'days_all_valid': len(a), 'mean_all_valid': round(a.pm.mean(), 1),
                      'geomean_all_valid': round(float(np.exp(np.log(a.pm).mean())), 1),
                      'days_with_blh': len(b), 'mean_days_with_blh': round(b.pm.mean(), 1),
                      'months_with_blh': ','.join(str(m) for m in sorted(b.d.dt.month.unique()))})
    # site availability inside the fixed panel
    s_ = sd[sd.city == city].assign(y=lambda q: q.d.dt.year)
    npanel = s_.station_id.nunique()
    for y, gy in s_[s_.y <= 2024].groupby('y'):
        per_day = gy.groupby('d').size()
        avail.append({'city': city, 'year': y, 'panel_stations': npanel, 'stations_with_any_day': gy.station_id.nunique(),
                      'days_with_data': len(per_day), 'mean_stations_per_day': round(per_day.mean(), 1),
                      'min_stations_per_day': int(per_day.min()), 'max_stations_per_day': int(per_day.max()),
                      'share_days_all_panel': round(float((per_day == npanel).mean()), 2)})

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, 'weather_adjustment_results.csv'), index=False)
pd.DataFrame(recon).to_csv(os.path.join(OUT, 'weather_adjustment_reconcile.csv'), index=False)
pd.DataFrame(avail).to_csv(os.path.join(OUT, 'weather_adjustment_availability.csv'), index=False)
json.dump({'seed': SEED, 'B': B, 'blh_gap': gap.astype(str).to_dict('records')},
          open(os.path.join(OUT, 'weather_adjustment.json'), 'w'), indent=1)
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
print(gap.to_string())
print(pd.DataFrame(recon).to_string())
print(res.to_string())
print(pd.DataFrame(avail).to_string())
