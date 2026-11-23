#!/usr/bin/env python3
"""Build the hourly and daily tables in data/raw/aq.duckdb from the XKDR Parquet files (run 01_download_xkdr.sh first).

Tables:
  st_state  station list with state and region (62 stations ship without a state; their names carry the city)
  pm        hourly PM2.5 and PM10, 0 < value <= 1500
  daily     station-day means, only days with at least 16 valid hours

Run:  python code/02_build_hourly_db.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import duckdb

# 62 stations ship without a state; their names carry the city, so map city -> state by hand.
CITY_STATE = {
    'Pune': 'Maharashtra', 'Navi Mumbai': 'Maharashtra',
    'Chennai': 'TamilNadu', 'Palkalaiperur': 'TamilNadu', 'Ooty': 'TamilNadu', 'Virudhunagar': 'TamilNadu',
    'Chengalpattu': 'TamilNadu', 'Karur': 'TamilNadu', 'Ariyalur': 'TamilNadu', 'Kanchipuram': 'TamilNadu',
    'Cuddalore': 'TamilNadu', 'Tirupur': 'TamilNadu', 'Dindigul': 'TamilNadu', 'Tirunelveli': 'TamilNadu',
    'Hosur': 'TamilNadu', 'Salem': 'TamilNadu', 'Tiruchirappalli': 'TamilNadu', 'Ranipet': 'TamilNadu',
    'Vellore': 'TamilNadu', 'Nagapattinam': 'TamilNadu',
    'Haveri': 'Karnataka', 'Hassan': 'Karnataka', 'Udupi': 'Karnataka', 'Hubballi': 'Karnataka',
    'Davanagere': 'Karnataka', 'Koppal': 'Karnataka', 'Raichur': 'Karnataka', 'Karwar': 'Karnataka',
    'Kalaburagi': 'Karnataka', 'Bidar': 'Karnataka', 'Gadag': 'Karnataka', 'Belgaum': 'Karnataka', 'Kolar': 'Karnataka',
    'Imphal': 'Manipur', 'Nandesari': 'Gujarat', 'Ahmedabad': 'Gujarat', 'Guwahati': 'Assam',
    'Kanpur': 'Uttar Pradesh', 'Baghpat': 'Uttar Pradesh', 'Lucknow': 'Uttar Pradesh', 'Shillong': 'Meghalaya',
    'Ernakulam': 'Kerala', 'Kozhikode': 'Kerala', 'Kochi': 'Kerala', 'Dhanbad': 'Jharkhand',
    'Pathardih': 'Jharkhand', 'Jorapokhar': 'Jharkhand', 'Agartala': 'Tripura', 'Ujjain': 'Madhya Pradesh',
    'Gwalior': 'Madhya Pradesh', 'Maihar': 'Madhya Pradesh', 'Indore': 'Madhya Pradesh', 'Suakati': 'Odisha',
    'Vijayawada': 'Andhra Pradesh', 'Ambala': 'Haryana', 'Darbhanga': 'Bihar',
}
REGION = {
    'North': ['Delhi', 'Haryana', 'Punjab', 'Chandigarh', 'Uttar Pradesh', 'Uttarakhand', 'Himachal Pradesh', 'Rajasthan'],
    'East': ['Bihar', 'Jharkhand', 'West Bengal', 'Odisha'],
    'Central': ['Madhya Pradesh', 'Chhattisgarh'],
    'West': ['Maharashtra', 'Gujarat'],
    'South': ['Karnataka', 'TamilNadu', 'Tamil Nadu', 'Kerala', 'Telangana', 'Andhra Pradesh', 'Puducherry'],
    'Northeast': ['Assam', 'Meghalaya', 'Sikkim', 'Arunachal Pradesh', 'Nagaland', 'Mizoram', 'Tripura', 'Manipur'],
}

c = duckdb.connect(str(HOURLY_DB))
c.sql(f"create or replace view m as select * from read_parquet('{RAW}/xkdr/measurements/**/*.parquet', hive_partitioning=true)")
c.sql(f"create or replace table st as select * from read_parquet('{RAW}/xkdr/stations.parquet')")

c.sql("create or replace table city_state(city varchar, state varchar)")
c.executemany("insert into city_state values (?, ?)", list(CITY_STATE.items()))
c.sql("create or replace table region_map(state varchar, region varchar)")
c.executemany("insert into region_map values (?, ?)", [(s, r) for r, ss in REGION.items() for s in ss])

# State: use the given one, else match a known city name inside the station name.
c.sql("""create or replace table st_state as
  select s.station_id, s.station_name, s.source, s.latitude lat, s.longitude lon,
    coalesce(nullif(s.state_name, 'NaN'),
             (select cs.state from city_state cs where s.station_name like '%' || cs.city || '%' limit 1)) state,
    coalesce(nullif(s.city_name, 'NaN'), nullif(s.city_name, '')) city
  from st s""")
c.sql("""update st_state set state = case
    when station_name like '%New Delhi%' then 'Delhi' when station_name like '%Mumbai%' then 'Maharashtra'
    when station_name like '%Chennai%' then 'TamilNadu' when station_name like '%Hyderabad%' then 'Telangana'
    when station_name like '%Kolkata%' then 'West Bengal' end
  where source = 'us_embassy'""")
c.sql("""create or replace table st_state as
  select s.*, coalesce(r.region, 'Other') region from st_state s left join region_map r on r.state = s.state""")

c.sql("""create or replace table pm as
  select m.station_id, m.source, s.state, s.region, coalesce(s.city, m.city_name) city,
         m.parameter_name p, m.collected_at ts, m.value v
  from m join st_state s using (station_id)
  where m.parameter_name in ('PM2.5', 'PM10') and m.value > 0 and m.value <= 1500""")

c.sql("""create or replace table daily as
  select station_id, source, state, region, city, p, cast(ts as date) d, avg(v) v, count(*) n
  from pm group by all having count(*) >= 16""")

print(c.sql("select region, count(distinct station_id) stations from st_state group by 1 order by 2 desc").fetchall())
print(c.sql("select count(*) from st_state where state is null").fetchall(), 'stations still without state')
print(c.sql("select p, count(*) station_days from daily group by 1").fetchall())
