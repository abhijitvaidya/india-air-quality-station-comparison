#!/usr/bin/env python3
"""Daily weather from the Open-Meteo historical archive (ERA5 reanalysis): rain, wind speed, temperature, relative humidity,
and the daily mean of hourly boundary-layer height, 2016-01-01 to 2026-10-08, for the cities of the weather adjustment.
Public, no key.

Run:  python code/10_download_weather.py
Output: data/raw/weather/city_day.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json
import os
import time
import urllib.request

OUT = os.path.dirname(WEATHER_CSV)
os.makedirs(OUT, exist_ok=True)
CITIES = {'Delhi': (28.61, 77.21), 'Varanasi': (25.32, 82.97), 'Agra': (27.18, 78.01), 'Kanpur': (26.45, 80.33),
          'Lucknow': (26.85, 80.95), 'Patna': (25.59, 85.14), 'Kolkata': (22.57, 88.36), 'Mumbai': (19.08, 72.88),
          'Bengaluru': (12.97, 77.59), 'Chennai': (13.08, 80.27), 'Hyderabad': (17.39, 78.49), 'Ludhiana': (30.90, 75.85),
          'Pune': (18.52, 73.86)}
out_p = os.path.join(OUT, 'city_day.csv'); done_p = os.path.join(OUT, 'done.txt')
done = set(open(done_p).read().split('\n')) if os.path.exists(done_p) else set()
if not os.path.exists(out_p):
    open(out_p, 'w').write('city,d,rain_mm,wind_kmh,temp_c,rh_pct,blh_m\n')
for city, (la, lo) in CITIES.items():
    if city in done: continue
    url = (f'https://archive-api.open-meteo.com/v1/archive?latitude={la}&longitude={lo}&start_date=2016-01-01&end_date=2026-10-08'
           '&daily=precipitation_sum,wind_speed_10m_mean,temperature_2m_mean,relative_humidity_2m_mean'
           '&hourly=boundary_layer_height&timezone=Asia%2FKolkata')
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as r: j = json.loads(r.read()); break
        except Exception as e:
            print(f'retry {city}: {e}', flush=True); time.sleep(60 * (attempt + 1))
    else:
        print(f'FAILED {city}', flush=True); continue
    blh = {}
    for t, v in zip(j['hourly']['time'], j['hourly']['boundary_layer_height']):
        if v is not None: blh.setdefault(t[:10], []).append(v)
    dly = j['daily']
    with open(out_p, 'a') as f:
        for i, d in enumerate(dly['time']):
            b = blh.get(d); bm = round(sum(b) / len(b)) if b else ''
            vals = [dly[k][i] for k in ('precipitation_sum', 'wind_speed_10m_mean', 'temperature_2m_mean', 'relative_humidity_2m_mean')]
            f.write(f'{city},{d},' + ','.join('' if v is None else str(v) for v in vals) + f',{bm}\n')
    open(done_p, 'a').write(city + '\n')
    print(f'{city} done ({len(dly["time"])} days)', flush=True)
    time.sleep(10)
print('ALL DONE', flush=True)
