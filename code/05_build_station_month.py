#!/usr/bin/env python3
"""Station-month means of continuous-station PM2.5 and PM10 from daily_all in data/raw/aq.duckdb (run 02 to 04 first).

A station-day is valid with at least 16 valid hours (daily table); each station-month keeps the mean of its valid days
and the count of valid days (ndays). The analyses use months with at least 15 valid days.
The repository already includes these files in data/derived/; this script rebuilds them from the raw downloads.

Run:  python code/05_build_station_month.py
Outputs: data/derived/continuous_station_month_pm25.csv, data/derived/continuous_station_month_pm10.csv
"""
import os, sys
import duckdb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *

c = duckdb.connect(str(HOURLY_DB), read_only=True)
for p, fn in [('PM2.5', 'continuous_station_month_pm25.csv'), ('PM10', 'continuous_station_month_pm10.csv')]:
    df = c.sql(f"""select station_id, strftime(d, '%Y%m') ym, avg(v) mean, count(*) ndays from daily_all
                  where p = '{p}' and source = 'cpcb_caaqm' and v > 0 group by all order by 1, 2""").df()
    df.to_csv(DERIVED / fn, index=False)
    print(fn, len(df), 'station-months')
