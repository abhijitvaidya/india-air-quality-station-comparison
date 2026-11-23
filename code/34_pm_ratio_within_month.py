#!/usr/bin/env python3
"""Standard deviation of the hourly ratio of PM2.5 to PM10 within each month: Ardhali Bazar by month, and stations across
India whose median monthly value is below 0.05 (Appendix G, Tables G4 and G5).

Run:  python code/34_pm_ratio_within_month.py   (needs data/raw/aq.duckdb; see data/README.md)
Outputs: results/pm_ratio_within_month.md, pm_ratio_monthly_ardhali.csv, pm_ratio_station_year.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os

import duckdb
import pandas as pd

HERE = str(RESULTS)
c = duckdb.connect(str(HOURLY_DB), read_only=True)
mo = c.sql("""
with h as (
  select station_id, date_trunc('hour', ts) ts,
         avg(case when p = 'PM2.5' then v end) a, avg(case when p = 'PM10' then v end) b
  from pm where source = 'cpcb_caaqm' and p in ('PM2.5', 'PM10') and ts >= timestamp '2018-01-01'
  group by 1, 2),
r as (select station_id, year(ts) y, month(ts) m, a / b r, ln(a) la, ln(b) lb from h where a > 0 and b > 0 and a / b <= 3)
select station_id, y, m, count(*) n, stddev_samp(r) sd, avg(r) mean_ratio, corr(la, lb) logcorr
from r group by 1, 2, 3 having count(*) >= 200
""").df()
names = c.sql("select station_id, station_name, state from st_state").df() if "st_state" in [t[0] for t in c.sql("show tables").fetchall()] else None
a = mo[mo.station_id == "site_273"].sort_values(["y", "m"])
a.round(4).to_csv(f"{HERE}/pm_ratio_monthly_ardhali.csv", index=False)
sy = mo.groupby(["station_id", "y"]).agg(months=("sd", "size"), median_sd=("sd", "median"), min_sd=("sd", "min"),
                                         median_logcorr=("logcorr", "median")).reset_index()
sy = sy[sy.months >= 6]
sy["fixed"] = sy.median_sd < 0.05
if names is not None:
    sy = sy.merge(names, on="station_id", how="left")
sy.round(4).to_csv(f"{HERE}/pm_ratio_station_year.csv", index=False)

L = ["# PM2.5/PM10 hourly ratio within months (flag A)", "",
     "Monthly SD of the hourly PM2.5/PM10 ratio; months with at least 200 paired hours. A station-year is 'fixed ratio' when the median monthly SD is below 0.05 (at least 6 months).", "",
     "## Ardhali Bazar, by month", "", "| year | " + " | ".join(str(i) for i in range(1, 13)) + " |", "|---|" + "---|" * 12]
for y, g in a.groupby("y"):
    d = g.set_index("m").sd
    L.append(f"| {y} | " + " | ".join(f"{d[i]:.3f}" if i in d.index else "" for i in range(1, 13)) + " |")
L += ["", "## Ardhali Bazar, by year", "", "| year | months | median monthly SD | lowest monthly SD | median log correlation |", "|---|---|---|---|---|"]
for _, r in sy[sy.station_id == "site_273"].iterrows():
    L.append(f"| {int(r.y)} | {int(r.months)} | {r.median_sd:.3f} | {r.min_sd:.3f} | {r.median_logcorr:.2f} |")
L += ["", "## All stations", "", "| year | station-years | median of median monthly SD | fixed ratio | share |", "|---|---|---|---|---|"]
for y, g in sy.groupby("y"):
    L.append(f"| {int(y)} | {len(g)} | {g.median_sd.median():.2f} | {int(g.fixed.sum())} | {100 * g.fixed.mean():.1f}% |")
L += ["", "## Fixed-ratio station-years", ""]
col = "station_name" if "station_name" in sy.columns else "station_id"
for y, g in sy[sy.fixed].groupby("y"):
    L.append(f"- {int(y)}: " + "; ".join(f"{str(n).split(' - ')[0]} ({s:.3f})" for n, s in zip(g[col], g.median_sd)))
open(f"{HERE}/pm_ratio_within_month.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
