#!/usr/bin/env python3
"""Satellite-derived monthly surface PM2.5 (WashU ACAG SatPM2.5 V6.GL.03, 0.01 degree, Asia files) at every station with
coordinates in data/derived/stations.csv, 2016-2024: the mean of the 5 x 5 cells around each station.
Source: https://www.satpm.org (public S3 bucket, no login). Each monthly file is about 92 MB; it is downloaded, read and
deleted. Resumable. The repository already includes the result in data/derived/.

Run:  python code/06_extract_satellite_pm25.py
Output: data/derived/satellite_pm25_station_month_v6.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import csv
import os
import subprocess

import netCDF4
import numpy as np

OUT = os.path.join(RAW, 'satpm'); os.makedirs(OUT, exist_ok=True)
CSV = os.path.join(DERIVED, 'satellite_pm25_station_month_v6.csv')
URL = 'https://satpmdata.s3.us-east-1.amazonaws.com/V6GL03/CNNPM25/Monthly/AS/V6GL03.CNNPM25.AS.{ym}-{ym}.nc'

with open(os.path.join(DERIVED, 'stations.csv')) as f:
    stations = [(r['station_id'], float(r['lat']), float(r['lon'])) for r in csv.DictReader(f) if r['lat']]

done = set()
if os.path.exists(CSV):
    with open(CSV) as f:
        done = {r['ym'] for r in csv.DictReader(f)}
else:
    with open(CSV, 'w') as f:
        f.write('station_id,ym,sat_pm25\n')

months = [f'{y}{m:02d}' for y in range(2016, 2025) for m in range(1, 13)]
for ym in months:
    if ym in done:
        continue
    path = os.path.join(OUT, f'{ym}.nc')
    r = subprocess.run(['curl', '-sf', '--retry', '3', '-m', '600', '-o', path, URL.format(ym=ym)])
    if r.returncode != 0:
        print(f'{ym} download failed ({r.returncode})', flush=True)
        continue
    d = netCDF4.Dataset(path)
    lat, lon, P = d.variables['lat'][:], d.variables['lon'][:], d.variables['PM25']
    rows = []
    for sid, la, lo in stations:
        i = int(np.argmin(np.abs(lat - la))); j = int(np.argmin(np.abs(lon - lo)))
        w = P[max(i - 2, 0):i + 3, max(j - 2, 0):j + 3]
        v = float(np.ma.mean(w)) if np.ma.count(w) else ''
        rows.append(f'{sid},{ym},{round(v, 2) if v != "" else ""}')
    d.close()
    os.remove(path)
    with open(CSV, 'a') as f:
        f.write('\n'.join(rows) + '\n')
    print(f'{ym} done ({len(rows)} stations)', flush=True)
print('ALL DONE', flush=True)
