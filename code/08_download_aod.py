#!/usr/bin/env python3
"""Monthly aerosol optical depth at 550 nm (CERES SYN1deg, through the NASA POWER API, parameter AOD_55), 2016-2024, at each
station location (rounded to 0.01 degree). POWER returns the value of the 1-degree cell that contains the point.
Public, no key. One request per location, paced.

Run:  python code/08_download_aod.py
Output: data/raw/power_aod55/aod55_<lat>_<lon>.json
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import datetime as dt, json, os, sys, time, urllib.request
import pandas as pd
from common import paired

CACHE = str(AOD_DIR); os.makedirs(CACHE, exist_ok=True)
URL = ('https://power.larc.nasa.gov/api/temporal/monthly/point?parameters=AOD_55&community=RE'
       '&longitude={lon}&latitude={lat}&start=2016&end=2024&format=JSON')
ids = paired('v6').station_id.unique()
st = pd.read_csv(f'{DERIVED}/stations.csv').dropna(subset=['lat', 'lon'])
st = st[st.station_id.isin(ids)]
pts = sorted({(round(a, 2), round(b, 2)) for a, b in zip(st.lat, st.lon)})
todo = [(la, lo) for la, lo in pts if not os.path.exists(f'{CACHE}/aod55_{la}_{lo}.json')]
print(len(pts), 'points;', len(todo), 'to fetch', flush=True)
manifest, fail = [], []
for n, (la, lo) in enumerate(todo, 1):
    f = f'{CACHE}/aod55_{la}_{lo}.json'; url = URL.format(lat=la, lon=lo); ok = False
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                b = r.read()
            j = json.loads(b); assert j['header'].get('sources') == ['SYN1DEG'], j['header'].get('sources')
            open(f, 'wb').write(b); ok = True; break
        except Exception as e:
            print('retry', la, lo, repr(e)[:120], flush=True); time.sleep(15)
    if ok:
        manifest.append({'file': os.path.basename(f), 'url': url, 'fetched_at_utc': dt.datetime.now(dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')})
    else:
        fail.append((la, lo))
    if n % 25 == 0: print(n, 'of', len(todo), flush=True)
    time.sleep(1.5)
if manifest:
    mp = f'{CACHE}/MANIFEST.csv'
    pd.DataFrame(manifest).to_csv(mp, mode='a', header=not os.path.exists(mp), index=False)
print('fetched', len(manifest), 'failed', fail)
