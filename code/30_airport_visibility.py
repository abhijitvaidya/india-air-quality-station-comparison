#!/usr/bin/env python3
"""Visibility at Varanasi, Patna, Lucknow and Delhi airports, 2016-2024: selected dry daytime reports, one observation per
clock hour (the median of its reports), and the shares below 2 km and 4 km (Appendix C, Table C1, Figure 4b).

Run:  python code/30_airport_visibility.py   (needs data/raw/metar; see data/README.md)
Outputs: results/airport_visibility_hourly.csv, results/airport_visibility_hourly.md
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import sys
import numpy as np, pandas as pd
from metar_codes import BAD, VIS

H = str(RESULTS)
rows = []
for icao in ['VEBN', 'VEPT', 'VILK', 'VIDP']:
    m = pd.read_csv(METAR_DIR / f'{icao}.csv', low_memory=False)
    m = m[(m.station == icao) & m.metar.notna()].copy()
    m['ist'] = pd.to_datetime(m.valid) + pd.Timedelta(hours=5, minutes=30)
    m['y'] = m.ist.dt.year; m = m[m.y.between(2016, 2024)]
    m['vis_m'] = m.metar.map(lambda s: (lambda g: int(g.group(1)) if g else np.nan)(VIS.search(' ' + s + ' ')))
    m['bad'] = m.metar.str.split(' RMK').str[0].map(lambda s: bool(BAD.search(s)))
    m['relh'] = pd.to_numeric(m.relh, errors='coerce')
    dup = int(m.duplicated(['valid']).sum())
    f = m[(m.vis_m > 0) & ~m.bad & (m.relh < 70) & m.ist.dt.hour.between(9, 17)].copy()
    f['t'] = f.ist.dt.floor('h'); f = f.sort_values('ist')
    for y, g in f.groupby('y'):
        hr = g.groupby('t').vis_m.median(); first = g.groupby('t').vis_m.first()
        slots = 9 * (366 if y % 4 == 0 else 365)
        rows.append({'icao': icao, 'year': y, 'reports': len(g), 'hours': len(hr), 'possible_hours': slots, 'reports_per_hour': round(len(g) / len(hr), 2),
                     'report_lt4000': (g.vis_m < 4000).mean(), 'hour_median_lt4000': (hr < 4000).mean(), 'hour_first_lt4000': (first < 4000).mean(),
                     'report_lt2000': (g.vis_m < 2000).mean(), 'hour_median_lt2000': (hr < 2000).mean(),
                     'report_ge5000': (g.vis_m >= 5000).mean(), 'hour_median_ge5000': (hr >= 5000).mean(),
                     'report_eq6000': (g.vis_m == 6000).mean(), 'hour_median_eq6000': (hr == 6000).mean(), 'hour_median_gt6000': (hr > 6000).mean(),
                     'duplicate_timestamps_all_reports': dup})
r = pd.DataFrame(rows); r.round(4).to_csv(f'{H}/airport_visibility_hourly.csv', index=False)
L = ['# Airport visibility with one observation per airport-hour', '']
for icao, name in [('VEBN', 'Varanasi'), ('VEPT', 'Patna'), ('VILK', 'Lucknow'), ('VIDP', 'Delhi')]:
    g = r[r.icao == icao]
    L += [f'## {name} ({icao})', '', '| Year | Selected reports | Distinct hours (of possible) | Below 4 km, per report | Below 4 km, per hour (median) | Below 4 km, per hour (first report) | Below 2 km, per hour | 5 km or more, per hour |', '|---|---|---|---|---|---|---|---|']
    for x in g.itertuples():
        L.append(f'| {x.year} | {x.reports:,} | {x.hours:,} ({x.possible_hours:,}) | {100 * x.report_lt4000:.0f}% | {100 * x.hour_median_lt4000:.0f}% | {100 * x.hour_first_lt4000:.0f}% | {100 * x.hour_median_lt2000:.0f}% | {100 * x.hour_median_ge5000:.0f}% |')
    L.append('')
open(f'{H}/airport_visibility_hourly.md', 'w').write('\n'.join(L)); print('\n'.join(L))
