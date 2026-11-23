"""Airport reports (METAR): station coordinates, weather codes that exclude a report, visibility parsing, and the
hourly median extinction for dry daytime hours (09-17 local time, relative humidity below 70%, no fog, mist,
precipitation, thunderstorm, dust or sand codes)."""
import os, re, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import METAR_DIR

AP = {'VEBN': ('Varanasi', 25.452, 82.859), 'VIDP': ('Delhi', 28.566, 77.103), 'VEPT': ('Patna', 25.591, 85.088),
      'VILK': ('Lucknow', 26.761, 80.889), 'VIAG': ('Agra', 27.157, 77.961), 'VECX': ('Kanpur', 26.404, 80.410),
      'VIAL': ('Prayagraj', 25.440, 81.734), 'VEGK': ('Gorakhpur', 26.740, 83.450), 'VABB': ('Mumbai', 19.089, 72.868),
      'VAAH': ('Ahmedabad', 23.077, 72.635), 'VECC': ('Kolkata', 22.654, 88.447),
      'VOBL': ('Bengaluru', 13.198, 77.706), 'VOHS': ('Hyderabad', 17.240, 78.429), 'VOMM': ('Chennai', 12.990, 80.169),
      'VIJP': ('Jaipur', 26.824, 75.812), 'VAPO': ('Pune', 18.582, 73.920), 'VILD': ('Ludhiana', 30.855, 75.952),
      'VIAR': ('Amritsar', 31.710, 74.797)}
BAD = re.compile(r'\b[-+]?(VC)?(MI|BC|PR|DR|BL|SH|TS|FZ)?(FG|BR|RA|DZ|SN|GR|GS|PL|DU|SA|SS|DS|PO|FC|UP)\b|\bTS\b')
VIS = re.compile(r'\s(\d{4})(NDV)?\s')


def load(icao):
    p = METAR_DIR / f'{icao}.csv'
    if p.exists():
        df = pd.read_csv(p, low_memory=False)
        return df[df.station == icao]
    return None


def dry_haze_hours(icao):
    """Hourly median extinction (bext, 1/km) for dry, hazy daytime hours at one airport; None if no METAR."""
    m = load(icao)
    if m is None or m.empty: return None
    m = m[m.metar.notna()].copy()
    m['vis_m'] = m.metar.map(lambda s: (lambda g: int(g.group(1)) if g else np.nan)(VIS.search(' ' + s + ' ')))
    m['bad'] = m.metar.map(lambda s: bool(BAD.search(s.split(' RMK')[0])))
    m['relh'] = pd.to_numeric(m.relh, errors='coerce')
    m['t'] = (pd.to_datetime(m.valid) + pd.Timedelta(hours=5, minutes=30)).dt.floor('h')
    m = m[(m.vis_m > 0) & ~m.bad & (m.relh < 70) & m.t.dt.hour.between(9, 17)]
    m['vis_km'] = np.minimum(m.vis_m, 10000) / 1000
    m['bext'] = 3.912 / m.vis_km
    return m.groupby('t').bext.median().reset_index()
