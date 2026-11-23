"""Shared loader: the paired table of continuous-station PM2.5 and satellite PM2.5 by station and month."""
import os
import sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DERIVED

DS = str(DERIVED)
SEASONS = ['winter', 'pre-monsoon', 'monsoon', 'post-monsoon']


def names():
    return pd.read_csv(f'{DS}/stations.csv')[['station_id', 'station_name', 'city', 'state', 'region']]


def cpcb_month(min_days=15, pol='pm25'):
    cp = pd.read_csv(f'{DS}/continuous_station_month_{pol}.csv', dtype={'ym': str})
    cp = cp[cp.ndays >= min_days].rename(columns={'mean': 'cpcb'})
    return cp[['station_id', 'ym', 'cpcb', 'ndays']]


def sat(version='v6'):
    return pd.read_csv(f'{DS}/satellite_pm25_station_month_{version}.csv', dtype={'ym': str}).dropna()


def paired(version='v6', min_days=15, years=(2016, 2024), cp=None):
    cp = cpcb_month(min_days) if cp is None else cp
    m = cp.merge(sat(version), on=['station_id', 'ym']).merge(names(), on='station_id')
    m = m[(m.cpcb > 0) & (m.sat_pm25 > 0)].copy()
    m['y'] = m.ym.str[:4].astype(int)
    m['mo'] = m.ym.str[4:].astype(int)
    if years:
        m = m[(m.y >= years[0]) & (m.y <= years[1])]
    m['season'] = np.select([m.mo.isin([11, 12, 1, 2]), m.mo.isin([3, 4, 5]), m.mo.isin([6, 7, 8, 9])],
                            ['winter', 'pre-monsoon', 'monsoon'], 'post-monsoon')
    m['lr'] = np.log(m.cpcb / m.sat_pm25)
    m['t'] = (m.y - 2016) * 12 + m.mo - 1
    return m.sort_values(['station_id', 't']).reset_index(drop=True)


def bh(p, q=0.05):
    """Benjamini-Hochberg rejections at false discovery rate q."""
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); thr = q * np.arange(1, n + 1) / n
    ok = np.where(p[o] <= thr)[0]
    out = np.zeros(n, bool)
    if len(ok): out[o[:ok.max() + 1]] = True
    return out


def by(p, q=0.05):
    """Benjamini-Yekutieli rejections (valid under any dependence)."""
    n = len(p)
    return bh(p, q / np.sum(1.0 / np.arange(1, n + 1)))
HERE = os.path.dirname(os.path.abspath(__file__))
