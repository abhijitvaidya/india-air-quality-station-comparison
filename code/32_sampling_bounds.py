#!/usr/bin/env python3
"""How high could the continuous-station mean be if the manual samples fell on polluted days? (Section 3.5, 4.6;
Tables D2 and D3)

Four Varanasi pairs of a manual station and the nearest continuous station. For each period, the continuous daily means
(valid day: at least 16 valid hours; a month is used with at least 15 valid days and a manual value) are
  - resampled 20,000 times on schedules of nine days a month,
  - reduced to the mean of the k highest days of each month (k = 9, 8, 6, 4, 2, 1).
The manual reports give no sample dates or counts, so the highest day of each month bounds any choice of days.

Needs the hourly database (data/raw/aq.duckdb; see data/README.md).
Run:  python code/32_sampling_bounds.py
Outputs: results/sampling_random_schedules.csv, results/sampling_bounds.csv, results/sampling_bhu_months.csv,
results/sampling_bounds.md
"""
import itertools, os, sys
import duckdb, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *

RNG = np.random.default_rng(20261115)
NDRAW = 20000
DAYS_PER_MONTH = 9          # 104 samples a year / 12
KS = [9, 8, 6, 4, 2, 1]
PAIRS = [("Banaras Hindu University", "site_5468", "IESD BHU", 0.6), ("Jawahar Nagar", "site_5469", "Bhelupur", 0.4),
         ("Saket Nagar", "site_5469", "Bhelupur", 1.5), ("Sigra", "site_5475", "Maldahiya", 1.3)]
WINDOWS = {"Nov 2023–Mar 2024": ("2023-11", "2024-03"), "Jul 2021–Mar 2022": ("2021-07", "2022-03"),
           "FY 2022–23": ("2022-04", "2023-03"), "FY 2023–24": ("2023-04", "2024-03"), "FY 2024–25": ("2024-04", "2025-03")}


def load_hourly(con, sid):
    h = con.execute("""select ts, v from pm where source='cpcb_caaqm' and p='PM10' and station_id=?
                       and v is not null and v > 0 order by ts""", [sid]).df()
    h["date"] = h.ts.dt.normalize(); h["hour"] = h.ts.dt.hour
    return h


def daily_table(h):
    g = h.groupby("date")
    d = pd.DataFrame({"n": g.v.size(), "mean24": g.v.mean()})
    d = d[d.n >= 16].copy()
    d["month"] = d.index.to_period("M"); d["dow"] = d.index.dayofweek
    return d


con = duckdb.connect(str(HOURLY_DB), read_only=True)
man = pd.read_csv(f'{DERIVED}/manual_station_month.csv')
man = man[(man.city == "Varanasi") & man.repeat_of.isna() & man.pm10.notna()].copy()
man["month"] = pd.PeriodIndex.from_fields(year=man.year, month=man.month, freq="M")
daily = {sid: daily_table(load_hourly(con, sid)) for sid in sorted({p[1] for p in PAIRS})}

win_rows, rows, mrows = [], [], []
for mst, sid, sname, km in PAIRS:
    d = daily[sid]; m_man = man[man.station == mst].set_index("month").pm10
    for wlabel, (s, e) in WINDOWS.items():
        months = [m for m in pd.period_range(s, e, freq="M") if m in m_man.index and (d.month == m).sum() >= 15]
        if len(months) < 3:
            continue
        dd = d[d.month.isin(months)]
        full = dd.groupby("month").mean24.mean()
        draws = np.zeros((NDRAW, len(months))); top9 = []
        for j, m in enumerate(months):
            x = dd[dd.month == m].mean24.to_numpy(); k = min(DAYS_PER_MONTH, len(x))
            idx = np.argsort(RNG.random((NDRAW, len(x))), axis=1)[:, :k]
            draws[:, j] = x[idx].mean(axis=1); top9.append(np.sort(x)[-k:].mean())
        rand_win = draws.mean(axis=1); man_win = m_man.reindex(months).mean(); cont_win = full.reindex(months).mean()
        win_rows.append({"manual_station": mst, "continuous_station": sname, "station_id": sid, "distance_km": km, "window": wlabel,
                         "months_used": len(months), "manual_mean": round(man_win, 1), "continuous_all_valid_days": round(cont_win, 1),
                         "random9_p2_5": round(np.percentile(rand_win, 2.5), 1), "random9_p97_5": round(np.percentile(rand_win, 97.5), 1),
                         "top9_days_each_month_bound": round(float(np.mean(top9)), 1),
                         "share_of_gap_random9_p97_5": round((np.percentile(rand_win, 97.5) - cont_win) / (man_win - cont_win), 2) if man_win > cont_win else np.nan,
                         "months_manual_above_top9": int(sum(m_man[m] > t for m, t in zip(months, top9)))})
        r = {'manual_station': mst, 'continuous_station': sname, 'window': wlabel, 'months': len(months), 'manual_mean': man_win,
             'continuous_mean': np.mean([d[d.month == m].mean24.mean() for m in months]),
             'valid_days_min': min(int((d.month == m).sum()) for m in months), 'valid_days_total': sum(int((d.month == m).sum()) for m in months)}
        for k in KS:
            r[f'top{k}'] = np.mean([np.sort(d[d.month == m].mean24.to_numpy())[-k:].mean() for m in months])
        rows.append(r)
        if mst == 'Banaras Hindu University' and wlabel.startswith('Nov 2023'):
            for m in months:
                x = np.sort(d[d.month == m].mean24.to_numpy())
                mrows.append({'month': str(m), 'manual': m_man[m], 'continuous_mean': x.mean(), 'valid_days': len(x), 'highest_day': x[-1], 'top9': x[-9:].mean(), 'top4': x[-4:].mean()})
wi = pd.DataFrame(win_rows); wi.to_csv(f'{RESULTS}/sampling_random_schedules.csv', index=False)
r = pd.DataFrame(rows); r.round(1).to_csv(f'{RESULTS}/sampling_bounds.csv', index=False)
b = pd.DataFrame(mrows); b.round(1).to_csv(f'{RESULTS}/sampling_bhu_months.csv', index=False)
L = ['# Sampling bounds, Varanasi pairs', '', '| Manual station | Continuous | Period | Manual mean | Continuous mean | Nine random days (95% of draws) | Valid days | ' +
     ' | '.join(f'k={k}' for k in KS) + ' |', '|' + '---|' * (7 + len(KS))]
for w, x in zip(wi.itertuples(), r.itertuples()):
    L.append(f'| {x.manual_station} | {x.continuous_station} | {x.window} | {x.manual_mean:.0f} | {x.continuous_mean:.0f} | {w.random9_p2_5:.0f} to {w.random9_p97_5:.0f} | '
             f'{x.valid_days_total} | ' + ' | '.join(f'{getattr(x, f"top{k}"):.0f}' for k in KS) + ' |')
L += ['', f'Station-months with the manual value above the nine-highest-day bound: {int(wi.months_manual_above_top9.sum())} of {int(wi.months_used.sum())}.',
      f'Highest day of each month below the manual mean in {int((r.top1 < r.manual_mean).sum())} of {len(r)} rows.',
      f'Banaras Hindu University, November 2023 to March 2024: manual {b.manual.mean():.1f}, continuous {b.continuous_mean.mean():.1f}; '
      f'November to February: {b.manual[:4].mean():.1f} and {b.continuous_mean[:4].mean():.1f}.']
open(f'{RESULTS}/sampling_bounds.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
