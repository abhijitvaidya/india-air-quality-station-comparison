"""File locations used by every script. Run scripts from the repository root, for example
    python code/11_national_series.py
"""
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / 'data' / 'derived'          # small derived inputs, included in the repository
RAW = ROOT / 'data' / 'raw'                  # downloads, not included (see data/README.md)
RESULTS = ROOT / 'results'
STEP = RESULTS / 'step_tests'
TABLES = ROOT / 'tables'
FIGURES = ROOT / 'figures'

CAMS_NC = RAW / 'cams' / 'eac4_monthly_india_2003-2025_pm.nc'
AOD_DIR = RAW / 'power_aod55'                 # one NASA POWER file per station location
METAR_DIR = RAW / 'metar'
WEATHER_CSV = RAW / 'weather' / 'city_day.csv'
HOURLY_DB = RAW / 'aq.duckdb'                 # hourly and daily database built from the XKDR files

for _d in (RESULTS, STEP, TABLES, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)


def aod_file(lat, lon):
    return AOD_DIR / f'aod55_{round(lat, 2)}_{round(lon, 2)}.json'
