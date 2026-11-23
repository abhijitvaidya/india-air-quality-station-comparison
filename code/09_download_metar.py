#!/usr/bin/env python3
"""Airport reports (METAR), 2016 to date, from the Iowa Environmental Mesonet ASOS/METAR archive (public, no key).
One request at a time, with pauses. Resumable per airport-year. The paper uses Varanasi (VEBN), Patna (VEPT),
Lucknow (VILK) and Delhi (VIDP).

Run:  python code/09_download_metar.py
Output: data/raw/metar/<ICAO>.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os
import time
import urllib.request

OUT = str(METAR_DIR); os.makedirs(OUT, exist_ok=True)
AIRPORTS = {  # ICAO: city
    'VEBN': 'Varanasi', 'VIAG': 'Agra', 'VILK': 'Lucknow', 'VECX': 'Kanpur', 'VIAL': 'Prayagraj', 'VEGK': 'Gorakhpur',
    'VIDP': 'Delhi', 'VEPT': 'Patna', 'VECC': 'Kolkata', 'VABB': 'Mumbai', 'VOBL': 'Bengaluru', 'VOHS': 'Hyderabad',
    'VOMM': 'Chennai', 'VAAH': 'Ahmedabad', 'VIJP': 'Jaipur', 'VAPO': 'Pune', 'VILD': 'Ludhiana', 'VIAR': 'Amritsar',
}
URL = ('https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station={s}&data=vsby&data=tmpf&data=dwpf&data=relh'
       '&data=sknt&data=wxcodes&data=metar&year1={y}&month1=1&day1=1&year2={y2}&month2=1&day2=1&tz=Etc/UTC'
       '&format=onlycomma&latlon=no&missing=M&trace=T&direct=no&report_type=3&report_type=4')
done_p = os.path.join(OUT, 'done.txt')
done = set(open(done_p).read().split()) if os.path.exists(done_p) else set()
for s, city in AIRPORTS.items():
    path = os.path.join(OUT, f'{s}.csv')
    for y in range(2016, 2027):
        key = f'{s}:{y}'
        if key in done: continue
        for attempt in range(4):
            try:
                with urllib.request.urlopen(URL.format(s=s, y=y, y2=y + 1), timeout=300) as r:
                    txt = r.read().decode('utf-8', 'replace'); break
            except Exception as e:
                print(f'retry {key}: {e}', flush=True); time.sleep(30 * (attempt + 1))
        else:
            print(f'FAILED {key}', flush=True); continue
        lines = txt.splitlines()
        if not lines: continue
        with open(path, 'a') as f:
            if os.path.getsize(path) == 0 if os.path.exists(path) else True:
                f.write(lines[0] + '\n')
            f.write('\n'.join(lines[1:]) + ('\n' if len(lines) > 1 else ''))
        open(done_p, 'a').write(key + '\n')
        print(f'{s} {city} {y}: {len(lines) - 1} reports', flush=True)
        time.sleep(5)
print('ALL DONE', flush=True)
