#!/usr/bin/env python3
"""Daily PM2.5 and PM10 for CPCB stations from OpenAQ, for the days the XKDR export lacks (January to March 2025 and from
2 September 2025). Only the 2025 PM10 values of the manual-station comparison (Table 4) use these days.

Needs an OpenAQ API key in the environment variable OPENAQ_API_KEY. Paced for OpenAQ's limit of 60 requests a minute.
Resumable: finished sensors are recorded in done.txt. With two dates, fetches only that window for the locations already
found, into data/raw/openaq/ext_<FROM>_<TO>/.

Run:  python code/03_download_openaq_2025.py            (2025-01-01 to 2026-10-03)
Outputs: data/raw/openaq/locations.parquet, sensor_*.parquet, days.parquet
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import json
import os
import sys
import time
import urllib.error
import urllib.request

import pandas as pd

BASE = os.path.join(RAW, 'openaq')
KEY = os.environ['OPENAQ_API_KEY']
EXTEND = len(sys.argv) == 3
FROM, TO = (sys.argv[1], sys.argv[2]) if EXTEND else ('2025-01-01', '2026-10-03')
OUT = os.path.join(BASE, f'ext_{FROM}_{TO}') if EXTEND else BASE
os.makedirs(OUT, exist_ok=True)


def get(path):
    for attempt in range(6):
        t0 = time.time()
        try:
            req = urllib.request.Request('https://api.openaq.org' + path, headers={'X-API-Key': KEY})
            with urllib.request.urlopen(req, timeout=60) as r:
                remaining = int(r.headers.get('X-Ratelimit-Remaining', 10))
                reset = int(r.headers.get('X-Ratelimit-Reset', 60))
                body = json.loads(r.read())
            dt = time.time() - t0
            if dt > 15:
                print(f'SLOW {dt:.0f}s {path[:90]}', flush=True)
            if remaining <= 1:
                print(f'PACE remaining={remaining} sleeping {reset + 1}s', flush=True)
                time.sleep(min(reset, 120) + 1)
            return body
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {'results': []}
            reset = int(e.headers.get('X-Ratelimit-Reset', 60) or 60)
            wait = min(reset, 300) + 2 if e.code == 429 else 5 * (attempt + 1)
            print(f'HTTP {e.code} after {time.time() - t0:.0f}s, reset={reset}, waiting {wait}s: {path[:90]}', flush=True)
            time.sleep(wait)
        except Exception as e:
            print(f'ERR {type(e).__name__} after {time.time() - t0:.0f}s: {str(e)[:80]} {path[:90]}', flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'failed: {path}')


locs = [l for l in get('/v3/locations?iso=IN&limit=1000')['results'] if (l.get('provider') or {}).get('name') == 'CPCB']
if EXTEND:  # same stations as the original backfill
    covered = set(pd.read_parquet(os.path.join(BASE, 'locations.parquet'))['location_id'])
    locs = [l for l in locs if l['id'] in covered]
pd.DataFrame([{'location_id': l['id'], 'name': l['name'], 'locality': l.get('locality'),
               'lat': l['coordinates']['latitude'], 'lon': l['coordinates']['longitude']} for l in locs]
             ).to_parquet(os.path.join(OUT, 'locations.parquet'))
sensors = [(l['id'], l['name'], s['id'], s['parameter']['name']) for l in locs for s in l['sensors']
           if s['parameter']['name'] in ('pm25', 'pm10')]

done_path = os.path.join(OUT, 'done.txt')
done = set(open(done_path).read().split()) if os.path.exists(done_path) else set()
print(f'{len(sensors)} sensors, {len(done)} already done', flush=True)

import threading
from concurrent.futures import ThreadPoolExecutor

# Shared pacing: at most ~50 requests/min across all workers (OpenAQ allows 60/min, 2000/hour).
_lock, _next = threading.Lock(), [0.0]
_get = get


def get(path):  # noqa: F811  (paced wrapper)
    with _lock:
        wait = _next[0] - time.time()
        if wait > 0:
            time.sleep(wait)
        _next[0] = max(time.time(), _next[0]) + 1.2
    return _get(path)


_done_lock = threading.Lock()


def fetch(item):
    i, (lid, name, sid, param) = item
    rows, page = [], 1
    while True:
        d = get(f'/v3/sensors/{sid}/days?date_from={FROM}&date_to={TO}&limit=1000&page={page}')
        res = d.get('results') or []
        for x in res:
            day = x['period']['datetimeFrom']['local'][:10]
            if FROM <= day <= TO:
                rows.append({'location_id': lid, 'sensor_id': sid, 'p': param, 'd': day, 'v': x['value'],
                             'n_hours': (x.get('coverage') or {}).get('observedCount')})
        if len(res) < 1000:
            break
        page += 1
        print(f'PAGE {page} sensor {sid} {name[:30]} {param}', flush=True)
    if rows:
        pd.DataFrame(rows).to_parquet(os.path.join(OUT, f'sensor_{sid}.parquet'))
    with _done_lock:
        with open(done_path, 'a') as f:
            f.write(f'{sid}\n')
    if i % 50 == 0:
        print(f'{i}/{len(sensors)} {name[:40]} {param} {len(rows)} days', flush=True)


todo = [(i, s) for i, s in enumerate(sensors) if str(s[2]) not in done]
with ThreadPoolExecutor(max_workers=4) as ex:
    list(ex.map(fetch, todo))

files = [os.path.join(OUT, f) for f in os.listdir(OUT) if f.startswith('sensor_')]
allr = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
allr.to_parquet(os.path.join(OUT, 'days.parquet'))
print(f'done: {len(allr)} sensor-days, {allr.location_id.nunique()} locations, {allr.d.min()} .. {allr.d.max()}', flush=True)
