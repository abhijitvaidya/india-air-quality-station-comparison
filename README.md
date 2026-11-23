# India air-quality stations compared with satellite, reanalysis and manual measurements

Code and derived data for the paper

> Vaidya, A. (2026). *Comparing India's air-quality stations with satellite, reanalysis and manual measurements,
> 2016–2024: larger gaps in Uttar Pradesh.* Preprint. https://doi.org/10.5281/zenodo.24445940

The paper and a copy of this repository (version 1.0.0) are archived together at Zenodo under that DOI.

The paper compares monthly PM2.5 at 543 of India's continuous monitoring stations with three references: satellite-derived
surface PM2.5 (WashU V6.GL.03), the CAMS EAC4 reanalysis and CERES SYN1deg aerosol optical depth. Each station is
compared with itself, calendar month by calendar month, across 2016–2024. It also compares PM10 at continuous stations
with CAMS, and PM10 at the manual stations of Uttar Pradesh from 2019 to 2025.

## What is here

| Folder | Contents |
|---|---|
| `code/` | Every script, numbered in the order they run |
| `data/derived/` | Small derived inputs: station list, station-month means, satellite PM2.5 at the stations, manual-station values. See `data/README.md` |
| `data/raw/` | Not included. Downloads land here; `data/README.md` says where each comes from |
| `results/` | The output of each analysis script, as CSV and markdown |
| `tables/`, `figures/` | The paper's tables and figures, made from `results/` |

## Reproducing the results

Python 3.14 with the packages in `requirements.txt`:

```
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

1. Download the references into `data/raw/` (`data/README.md`): CAMS by hand (`code/07_download_cams.md`), then
   `python code/08_download_aod.py`, `python code/09_download_metar.py` and `python code/10_download_weather.py`.
2. The checks at Ardhali Bazar (scripts 29 and 31 to 34) also need the hourly database. Build it with
   `code/01_download_xkdr.sh`, `code/02_build_hourly_db.py`, `code/03_download_openaq_2025.py` and `code/04_merge_openaq.py`.
   This needs API keys from XKDR and OpenAQ, set as the environment variables `XKDR_KEY` and `OPENAQ_API_KEY`.
3. Run `bash run_all.sh`. It runs scripts 11 to 41 and rewrites `results/`, `tables/` and `figures/`.

The PM2.5 and PM10 comparisons (scripts 11 to 22) need only `data/derived/`, the CAMS file and the optical-depth files.
Scripts 05 and 06 rebuild `data/derived/` from the raw downloads; the included files are what they produced.

## Scripts

| Script | What it does | Paper |
|---|---|---|
| `01_download_xkdr.sh`, `02_build_hourly_db.py` | Hourly CPCB data from the XKDR Forum database; hourly and daily tables | Section 2 |
| `03_download_openaq_2025.py`, `04_merge_openaq.py` | Days missing from the XKDR export in 2025, from OpenAQ | Section 2 |
| `05_build_station_month.py` | Station-month means | Section 2 |
| `06_extract_satellite_pm25.py` | Satellite PM2.5 at each station | Section 2 |
| `07_download_cams.md`, `08` to `10` | CAMS, optical depth, airport reports, weather | Section 2 |
| `11_national_series.py` | Each station against itself, year by year, against three references | 3.1, 4.1, 4.2; Table 2, E1 |
| `12_national_monthly.py` | Monthly series, linked check, steady-ratio test | Figure 1, Table 3, E3 |
| `13_single_station.py` | Two-period change at each eligible station | 3.3, 4.4; Table F1, Figure 2 |
| `14_gap_by_year.py` | Uttar Pradesh minus the rest of India, by year | 4.2; Table E2 |
| `15_pm25_base_periods.py`, `16_single_base_years.py` | Other base periods | 4.2; Table A3 |
| `17_single_station_checks.py`, `18_period_choice.py` | Eligibility rules, years left out, period choice, Ardhali Bazar by year | Appendix A2; Tables A1, A2; Figure 3 |
| `19_pm10_manual.py`, `20_pm10_manual_checks.py` | PM10 against CAMS; manual and continuous stations | 3.5, 4.6; Table 4 |
| `21_pm10_base_periods.py`, `22_pm25_pm10_ratio.py` | PM10 base periods; PM2.5 and PM10 at the same stations | 4.2; Tables A3, A4 |
| `23_step_tests.py` to `25_trend_and_dependence.py` | Step and trend models, dependence between stations | 3.4; Tables A6 to A8 |
| `26_interval_coverage.py`, `27_achieved_coverage.py` | Coverage of the single-station intervals and the widening factor | Appendix A5; Table A5 |
| `29_weather_adjustment.py` | Weather adjustment | 3.6, 4.5; Table B1 |
| `30_airport_visibility.py`, `31_airport_paired_pm25.py` | Airport visibility | 4.5; Table C1 |
| `32_sampling_bounds.py` | Manual and continuous PM10 under different sampling days | 3.5, 4.6; Tables D2, D3 |
| `33_ardhali_record_checks.py`, `34_pm_ratio_within_month.py` | Checks on the Ardhali Bazar record | 4.5; Appendix G |
| `40_make_tables.py`, `41_make_figures.py` | Tables and figures | |

Shared modules: `paths.py` (file locations), `common.py` (the paired station and satellite table), `two_period.py`
(two-period change at one station), `metar_codes.py` (airport report filters).

## Licence

Code: MIT (`LICENSE`). Data: see `DATA_LICENCE.md`. The manual-station values come from UPPCB's reports and remain under
UPPCB's terms.

## Citation

See `CITATION.cff`.
