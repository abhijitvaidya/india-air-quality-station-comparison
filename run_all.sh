#!/bin/bash
# Reproduce every result, table and figure from data/derived/ and the downloads in data/raw/.
# Run from the repository root:  bash run_all.sh
# Steps 01 to 10 rebuild data/derived/ and fetch data/raw/ (see data/README.md); they are not run here.
set -euo pipefail
cd "$(dirname "$0")"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$(mktemp -d)}"
PY="${PYTHON:-python3}"
run() { echo "== $*"; "$PY" -W ignore "$@" > /dev/null; }

# PM2.5 against the three references (needs data/raw/cams and data/raw/power_aod55)
run code/11_national_series.py
run code/12_national_monthly.py
run code/13_single_station.py
run code/14_gap_by_year.py
run code/15_pm25_base_periods.py
run code/16_single_base_years.py
run code/17_single_station_checks.py
run code/18_period_choice.py
# PM10 and the manual stations
run code/19_pm10_manual.py
run code/20_pm10_manual_checks.py
run code/21_pm10_base_periods.py
run code/22_pm25_pm10_ratio.py
# Step and trend models
for v in orig month12 calendar sync trim24 drop_2020 covid_dummy combined trend combined_trend resid_null; do
  run code/23_step_tests.py "$v"
done
run code/23_step_tests.py orig v5
run code/23_step_tests.py trend v5
run code/24_step_summary.py
run code/25_trend_and_dependence.py
# Interval coverage (about an hour in all)
for k in 0 1 2; do run code/26_interval_coverage.py "$k"; done
run code/26_interval_coverage.py
run code/27_achieved_coverage.py
# Checks at Ardhali Bazar (need data/raw/aq.duckdb, data/raw/weather and data/raw/metar)
run code/29_weather_adjustment.py
run code/30_airport_visibility.py
run code/31_airport_paired_pm25.py
run code/32_sampling_bounds.py
run code/33_ardhali_record_checks.py
run code/34_pm_ratio_within_month.py
# Tables and figures
run code/40_make_tables.py
run code/41_make_figures.py
echo "done"
