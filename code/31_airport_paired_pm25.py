#!/usr/bin/env python3
"""Varanasi airport: selected dry daytime hours paired with the median PM2.5 of continuous stations within 30 km that
report at least 300 such hours in both 2019 and 2024 (Ardhali Bazar is the only one), by year (Section 4.5, Table C1,
Figure 4c).

Run:  python code/31_airport_paired_pm25.py   (needs data/raw/metar and data/raw/aq.duckdb; see data/README.md)
Output: results/airport_paired_pm25.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, re, shutil, sys
import duckdb, numpy as np, pandas as pd
import metar_codes as mc
city, la, lo = mc.AP['VEBN']
c = duckdb.connect(str(HOURLY_DB), read_only=True)
FIXED_SQL = """and pm.station_id in (select station_id from pm where p='PM2.5' and source='cpcb_caaqm' and hour(ts) between 9 and 17
    group by 1 having count(*) filter (where year(ts)=2019) >= 300 and count(*) filter (where year(ts)=2024) >= 300)"""
cp = c.sql(f"""select date_trunc('hour', ts) t, median(v) pm, count(*) n_st, string_agg(distinct station_id) ids from pm join st_state s using (station_id)
  where pm.p='PM2.5' and pm.source='cpcb_caaqm' and s.lat is not null and hour(ts) between 9 and 17
    and 2*6371*asin(sqrt(pow(sin(radians(s.lat-{la})/2),2)+cos(radians({la}))*cos(radians(s.lat))*pow(sin(radians(s.lon-{lo})/2),2))) <= 30
    {FIXED_SQL} group by 1 order by 1""").df()
h = mc.dry_haze_hours('VEBN')
j = h.merge(cp[['t', 'pm']], on='t'); j['y'] = j.t.dt.year
rows = []
for y, g in j.groupby('y'):
    if len(g) < 150: continue
    vis = 3.912 / g.bext
    rows.append({'y': y, 'hours': len(g), 'vis_km_med': round(3.912 / g.bext.median(), 1), 'cpcb_pm_med': round(g.pm.median()),
                 'cpcb_pm_med_exact': g.pm.median(), 'share_below_4km': (vis < 4).mean(), 'corr_log': round(np.corrcoef(np.log(g.bext), np.log(g.pm))[0, 1], 2)})
r = pd.DataFrame(rows); r.to_csv(f'{RESULTS}/airport_paired_pm25.csv', index=False)
print('fixed stations:', sorted(set(','.join(cp.ids).split(','))))
print(r.to_string(index=False))
