#!/usr/bin/env python3
"""Build daily_all in data/raw/aq.duckdb: XKDR daily values, plus OpenAQ daily values for the days XKDR lacks.
OpenAQ locations are matched to XKDR stations by coordinates (within 300 m). In the April to August 2025 overlap the
median OpenAQ/XKDR ratio is 0.999 (PM2.5) and 0.997 (PM10). Two stations that disagree are left out of the fill:
Lodhi Road (IMD), ratio 0.79, and Jodhpur Collectorate, correlation 0.51. Where a sensor-day appears in more than one
fetch, the latest fetch wins.

Run after 02 and 03:  python code/04_merge_openaq.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import glob
import duckdb

EXCLUDE = ('Lodhi Road, New Delhi - IMD', 'Collectorate, Jodhpur - RSPCB')

c = duckdb.connect(str(HOURLY_DB))
# (file, batch): batch 0 = original fetch, then ext_ dirs in name (date) order.
files = [(f, 0) for f in glob.glob(f'{RAW}/openaq/sensor_*.parquet')]
for i, ext in enumerate(sorted(glob.glob(f'{RAW}/openaq/ext_*/')), start=1):
    files += [(f, i) for f in glob.glob(ext + 'sensor_*.parquet')]
batch = ' union all '.join(f"select *, {b} batch from read_parquet({[f for f, bb in files if bb == b]})"
                           for b in sorted({b for _, b in files}))
c.sql(f"""create or replace table oq as
  select * exclude (batch) from ({batch})
  qualify row_number() over (partition by sensor_id, d order by batch desc) = 1""")
c.sql(f"create or replace table oloc as select * from read_parquet('{RAW}/openaq/locations.parquet')")
c.sql("""create or replace table oq_map as
  select o.location_id, o.name oname, s.station_id,
    2*6371*asin(sqrt(pow(sin(radians(s.lat-o.lat)/2),2)+cos(radians(o.lat))*cos(radians(s.lat))*pow(sin(radians(s.lon-o.lon)/2),2))) km
  from oloc o, st_state s where s.source='cpcb_caaqm' and s.lat is not null
  qualify row_number() over (partition by o.location_id order by km)=1""")

c.sql(f"""create or replace table daily_all as
  select station_id, source, state, region, city, p, d, v, n, 'xkdr' src from daily
  union all
  select s.station_id, s.source, s.state, s.region, s.city,
         case o.p when 'pm25' then 'PM2.5' else 'PM10' end p, cast(o.d as date) d, avg(o.v) v,
         max(o.n_hours) n, 'openaq' src
  from oq o join oq_map m on m.location_id = o.location_id and m.km <= 0.3 and m.oname not in {EXCLUDE}
  join st_state s on s.station_id = m.station_id
  where o.n_hours >= 16
    and not exists (select 1 from daily x where x.station_id = s.station_id and x.d = cast(o.d as date)
                    and x.p = case o.p when 'pm25' then 'PM2.5' else 'PM10' end)
  group by all""")

print(c.sql("""select strftime(d,'%Y-%m') ym, count(distinct station_id) n_stations,
    count(*) filter (where src='xkdr') xkdr_days, count(*) filter (where src='openaq') openaq_days
  from daily_all where p='PM2.5' and d >= '2024-10-01' group by 1 order by 1""").df().to_string(index=False))
