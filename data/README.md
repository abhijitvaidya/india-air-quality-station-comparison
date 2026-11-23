# Data

## Included: `data/derived/`

| File | Rows | Columns | Made by |
|---|---|---|---|
| `stations.csv` | 553 | station_id, station_name, city, state, region, lat, lon, coord_source | XKDR station list; coordinates for 59 stations added from OpenAQ (40) and the CPCB station listing (19); 3 stations have none |
| `continuous_station_month_pm25.csv` | 35,226 | station_id, ym, mean (µg/m³), ndays (valid days) | `code/05_build_station_month.py` |
| `continuous_station_month_pm10.csv` | 33,672 | as above | `code/05_build_station_month.py` |
| `satellite_pm25_station_month_v6.csv` | 59,940 | station_id, ym, sat_pm25 (µg/m³; mean of 5 × 5 cells) | `code/06_extract_satellite_pm25.py` |
| `satellite_pm25_station_month_v5.csv` | 53,568 | as above, earlier version V5.GL.06 | the same extraction on V5.GL.06 files |
| `manual_station_month.csv` | 6,324 | year, month, city, station, pm10, so2, no2, pm25 (µg/m³), repeat_of | UPPCB monthly reports, January 2019 to December 2025 |
| `manual_city_centres.csv` | 19 | city, lat, lon, source | OpenStreetMap Nominatim, query "<city>, Uttar Pradesh, India" |
| `uppcb_reports.csv` | 84 | month, url, sha256, retrieved | the report addresses |

Notes:

- A station-day is valid with at least 16 valid hours. The analyses use months with at least 15 valid days.
- `coord_source` is `database` for coordinates from the XKDR station list, and `openaq` or `cpcb_portal` for the 59
  stations whose coordinates were added. The interval coverage simulation uses the eligible stations with `database`
  coordinates.
- Continuous data to 1 September 2025 come from the XKDR export of 30 September 2026. Days from January to March 2025
  and from 2 September 2025 come from OpenAQ's copy of the same CPCB feed. They enter only the 2025 PM10 values of the
  manual-station comparison.
- `repeat_of` marks a station-month whose PM10 and SO2 both repeat an earlier month at the same station (92 rows).
  These are left out as possible reprints. PM2.5 appears in the reports from July 2024 only.
- The manual reports give no station coordinates. A manual station is compared with CAMS at the mean location of its
  city's continuous stations, or at the city centre for the 19 cities without one.

## Not included: `data/raw/`

| Path | Source | Access | Script |
|---|---|---|---|
| `xkdr/` | XKDR Forum India air quality database, https://airquality.xkdr.org | free API key (`XKDR_KEY`) | `code/01_download_xkdr.sh` |
| `aq.duckdb` | built from `xkdr/` and `openaq/` | | `code/02_build_hourly_db.py`, `code/04_merge_openaq.py` |
| `openaq/` | OpenAQ API v3, https://openaq.org | free API key (`OPENAQ_API_KEY`) | `code/03_download_openaq_2025.py` |
| `satpm/` (temporary) | WashU SatPM2.5 V6.GL.03 monthly files, https://www.satpm.org | public | `code/06_extract_satellite_pm25.py` |
| `cams/eac4_monthly_india_2003-2025_pm.nc` | CAMS EAC4 monthly averaged fields, Atmosphere Data Store | free account | `code/07_download_cams.md` |
| `power_aod55/` | CERES SYN1deg aerosol optical depth through NASA POWER, https://power.larc.nasa.gov | public | `code/08_download_aod.py` |
| `metar/` | Iowa Environmental Mesonet ASOS/METAR archive, https://mesonet.agron.iastate.edu | public | `code/09_download_metar.py` |
| `weather/city_day.csv` | Open-Meteo historical archive (ERA5), https://open-meteo.com | public | `code/10_download_weather.py` |

Versions used in the paper:

| Source | Version | Retrieved |
|---|---|---|
| XKDR Forum database | export of 30 September 2026 | October 2026 |
| OpenAQ | API v3 daily values | October 2026 |
| WashU | V6.GL.03 (V5.GL.06 for one check) | October and November 2026 |
| CAMS EAC4 | monthly means 2003–2025; SHA-256 6a3b62f815f2c3dd87c7673dad17e38f8328b7793148cbac4f6c2a5abacd984d | 17 November 2026 |
| CERES SYN1deg | Edition 4.2, NASA POWER API 2.10.0 | 18 to 21 November 2026 |
| METAR | Iowa Environmental Mesonet archive | October 2026 |
| Weather | Open-Meteo historical archive | October 2026 |

Later downloads can differ: XKDR and OpenAQ may revise past values, and CAMS and NASA POWER may reprocess. The
included `data/derived/` files and `results/` are the versions behind the paper.
