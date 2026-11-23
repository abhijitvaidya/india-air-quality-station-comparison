#!/usr/bin/env python3
"""Tables of the paper as markdown fragments, from the result files (Table 2 and the appendix tables; Tables 1, 3, 4
and A2 to A4 are typed in the paper and their numbers are in the result files named in tables/README.md).

Run after the analysis scripts:  python code/40_make_tables.py
Writes tables/*.md. Reads only; changes no result file.
"""
import hashlib, math, os, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *

OUT = TABLES
for old in OUT.glob("*.md"):
    if old.name != "README.md":
        old.unlink()


def m(x, d=0, plus=False):
    """Number with a true minus sign."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    if round(x, d) == 0:
        x = 0.0   # no "-0"
    s = f"{x:+.{d}f}" if plus else f"{x:.{d}f}"
    return s.replace("-", "−")


pc = lambda x: m(x, 0, plus=True) + "%"
iv = lambda lo, hi, unit="%": f"({m(lo, 0, plus=True)}{unit} to {m(hi, 0, plus=True)}{unit})"


def table(name, head, rows):
    L = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    L += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    (OUT / f"{name}.md").write_text("\n".join(L) + "\n")
    print(name, len(rows), "rows")


def short(n):
    return re.sub(r"\s+-\s+[^-]*$", "", str(n)).strip()   # drop only the trailing agency name


def md_table(path, header):
    """Rows of the first markdown table after a heading that starts with `header`."""
    L = path.read_text().splitlines()
    i = next(k for k, x in enumerate(L) if x.startswith(header))
    j = i + 1
    while not L[j].startswith("|"): j += 1
    rows = []
    while j < len(L) and L[j].startswith("|"):
        rows.append([c.strip() for c in L[j].strip("|").split("|")]); j += 1
    return rows[0], rows[2:]


REF = {"WashU": "Satellite PM2.5", "CAMS": "CAMS PM2.5", "CERES": "Aerosol optical depth"}
REFK = {"washu": "Satellite PM2.5", "cams": "CAMS PM2.5", "ceres": "Aerosol optical depth"}

# ------------------------------------------------------------------ Table 2 (main text): groups against each reference
_, body = md_table(RESULTS / "national_series.md", "## Mean of 2023 and 2024")
keep = ["Uttar Pradesh", "Rest of India", "Plain outside Uttar Pradesh"]
rows = [[r[0], REF[r[1]], r[2], r[3], r[4], r[5]] for g in keep for r in body if r[0] == g]
assert len(rows) == 9
table("T2_groups", ["Group", "Reference", "Stations", "Station PM2.5", "Reference at the stations", "Station against reference (95% interval)"], rows)

# ------------------------------------------------------------------ A1 two-period change under other eligibility rules (166 stations)
e = pd.read_csv(RESULTS / "two_period_variants.csv").set_index("variant")
lab = {"orig": "Main rule", "drop_partial_winters": "Winters cut by a period boundary dropped", "winter_calendar_year": "Winter assigned by calendar year",
       "all_four_seasons": "All four seasons required", "min3_per_season": "At least three months a season", "common_calendar_months": "Same calendar months in both periods"}
rows = [[lab[k], int(e.stations[k]), int(e.up[k]), m(e.median_change[k], 1) + "%", m(e.ardhali_change[k], 1) + "%", int(e.ardhali_rank[k])] for k in lab]
table("A1_rules", ["Rule", "Stations", "In Uttar Pradesh", "Median change", "Ardhali Bazar", "Ardhali Bazar rank"], rows)

# ------------------------------------------------------------------ A5 coverage and widening (simulation on 159 of the 166 eligible stations)
cov = pd.read_csv(RESULTS / "interval_coverage.csv").set_index("scenario")
wid = {k: [k, "", f"{cov.widening_factor[k]:.2f}"] for k in cov.index}
ach = pd.read_csv(RESULTS / "achieved_coverage.csv").set_index("scenario")   # fresh draws, factor 1.62 throughout
lab = {"no year effect": "No year-to-year variation", "estimated year share": "Year-to-year variation as estimated (15% of variance)",
       "double year share": "Twice the estimated year-to-year variation"}
rows = [[lab[k], f"{100 * cov.loc[k, 'coverage_block_corr']:.0f}%", wid[k][2], f"{100 * ach.loc[k, 'widened']:.1f}%"] for k in lab]
assert wid["estimated year share"][2] == "1.62"
table("A5_coverage", ["Assumption in the simulation", "Coverage of the nominal 95% interval", "Widening factor", "Achieved coverage with the factor of 1.62"], rows)

# ------------------------------------------------------------------ A6 step variants, A7 best model, A8 residual correlation (175 stations)
v = pd.read_csv(STEP / "step_test_variants.csv").set_index("variant")
names = [("orig", "Single step, four season terms (main specification)"), ("month12", "Twelve month terms in place of four seasons"),
         ("calendar", "Bootstrap blocks of six calendar months"), ("sync", "Same block dates at every station"),
         ("trim24", "Step at least 24 months from either end"), ("drop_2020", "2020 excluded"),
         ("covid_dummy", "Separate term for March to June 2020"),
         ("combined", "Blocks of 12 calendar months, 12-month trim, March to June 2020 excluded"),
         ("trend", "Linear trend allowed (step beyond the trend)"), ("combined_trend", "Combined specification with linear trend"),
         ("resid_null", "Main specification, residuals drawn from the fit without a step")]
rows = []
for k, nm in names:
    r = v.loc[k]
    assert int(r.tested) == 175
    rows.append([nm, int(r.bh_5pct), int(r.by_5pct), f"{int(r.up_down)} / {int(r.up_up)}", f"{int(r.out_down)} / {int(r.out_up)}",
                 f"{m(r.ardhali_step_pct)}%, p = {r.ardhali_p:.3f}"])
table("A6_step_variants", ["Specification", "Detected, BH", "Detected, BY", "Uttar Pradesh, down / up", "Elsewhere, down / up", "Ardhali Bazar step"], rows)
t = pd.read_csv(STEP / "trend_vs_step.csv")
rows = []
for nm, g in [("Uttar Pradesh", t[t.state == "Uttar Pradesh"]), ("Elsewhere", t[t.state != "Uttar Pradesh"])]:
    c = g.best_bic.value_counts()
    rows.append([nm, len(g), int(c.get("season only", 0)), int(c.get("trend", 0)), int(c.get("step", 0)), int(c.get("trend+step", 0)),
                 f"{int((g.trend_pct_per_year < 0).sum())}", m(g.trend_pct_per_year.median(), 1)])
table("A7_best_model", ["Group", "Stations", "Seasons only", "Trend", "Step", "Trend and step", "Negative trend", "Median trend, % a year"], rows)
print("A7 check (main text: 18 of 21, −4.6; 84 of 154, −0.5):", rows[0][1], rows[0][6], rows[0][7], rows[1][1], rows[1][6], rows[1][7])
d = pd.read_csv(STEP / "residual_correlation_by_distance.csv")
table("A8_resid_corr", ["Distance between stations, km", "Station pairs", "Median correlation"],
      [[b.replace("-", "–").replace("600+", "over 600"), int(n), f"{md:.2f}"] for b, n, md in zip(d.band, d["size"], d["median"])])

# ------------------------------------------------------------------ B1 weather
w = pd.read_csv(RESULTS / "weather_adjustment_results.csv")
av = pd.read_csv(RESULTS / "weather_adjustment_availability.csv").groupby("city").panel_stations.max()
C = w[w.spec.str.startswith("C ")].set_index("city"); E = w[w.spec.str.startswith("E ")].set_index("city")
rows = []
for city in ["Varanasi", "Delhi", "Patna", "Lucknow", "Kanpur"]:
    r = C.loc[city]
    rows.append([city + (" (Ardhali Bazar)" if city == "Varanasi" else ""), int(av[city]), m(r.raw_geo_2023) + "%",
                 f"{m(r.weather_2023)}% ({m(r.weather_2023_lo)} to {m(r.weather_2023_hi)})", m(r.raw_geo_2024) + "%",
                 f"{m(r.weather_2024)}% ({m(r.weather_2024_lo)} to {m(r.weather_2024_hi)})", m(E.loc[city, "weather_2024"]) + "%"])
table("B1_weather", ["City", "Stations in panel", "2019 to 2023, unadjusted", "2019 to 2023, adjusted", "2019 to 2024, unadjusted", "2019 to 2024, adjusted",
                     "2019 to 2024, adjusted, station fixed effects"], rows)

# ------------------------------------------------------------------ C1 airport (one observation per airport-hour)
ah = pd.read_csv(RESULTS / "airport_visibility_hourly.csv")
pr = pd.read_csv(RESULTS / "airport_paired_pm25.csv").set_index("y")
vb = ah[ah.icao == "VEBN"].set_index("year"); pt = ah[ah.icao == "VEPT"].set_index("year")
rows = []
for y in range(2016, 2025):
    assert pr.hours[y] <= vb.hours[y] <= vb.possible_hours[y]
    rows.append([y, f"{int(vb.reports[y]):,}", f"{int(vb.hours[y]):,}", f"{100 * vb.hour_median_lt2000[y]:.0f}%", f"{100 * vb.hour_median_lt4000[y]:.0f}%",
                 f"{100 * vb.hour_median_ge5000[y]:.0f}%", f"{100 * vb.hour_median_eq6000[y]:.0f}%", f"{int(pr.hours[y]):,}", int(pr.cpcb_pm_med[y]),
                 f"{100 * pt.hour_median_lt4000[y]:.0f}%"])
table("C1_airport", ["Year", "Selected reports", "Selected hours", "Below 2 km", "Below 4 km", "5 km or more", "At 6,000 m", "Hours paired with the station",
                     "Station PM2.5 in paired hours, median", "Patna airport, below 4 km"], rows)

# ------------------------------------------------------------------ D1 manual stations by city
man = pd.read_csv(DERIVED / "manual_station_month.csv")
man = man[man.pm10.notna() & (man.pm10 > 0)]
nrep = int(man.repeat_of.notna().sum())
man = man[man.repeat_of.isna()]
st = pd.read_csv(DERIVED / "stations.csv")
p10 = pd.read_csv(DERIVED / "continuous_station_month_pm10.csv", dtype={"ym": str})
p10 = p10[(p10.ndays >= 15) & (p10["mean"] > 0) & (p10.ym >= "201901") & (p10.ym <= "202512")]
upc = st[(st.state == "Uttar Pradesh") & st.lat.notna() & st.station_id.isin(p10.station_id)]
alias = {"Prayagraj": ["Prayagraj", "Allahabad"], "Raebareli": ["Raebareli", "Rae Bareli"]}
rows = []; both = 0
for city, g in man.groupby("city"):
    nc = int(upc.city.isin(alias.get(city, [city])).sum()); both += nc > 0
    rows.append([city, g.station.nunique(), (f"{int(g.year.min())}–{int(g.year.max())}" if g.year.min() != g.year.max() else str(int(g.year.min()))), f"{len(g):,}", nc if nc else "none"])
rows.append(["All", man.groupby(["city", "station"]).ngroups, "2019–2025", f"{len(man):,}", int(sum(r[4] for r in rows if r[4] != "none"))])
table("D1_manual_cities", ["City", "Manual stations", "Years with data", "Station-months of PM10", "Continuous stations with PM10"], rows)
print("D1 check: cities", man.city.nunique(), "stations", rows[-1][1], "values", len(man), "reprints left out", nrep, "| cities with both networks", both, "continuous stations", rows[-1][4])

# ------------------------------------------------------------------ D2 sampling bound, D3 the university pair by month
s = pd.read_csv(RESULTS / "sampling_random_schedules.csv")
fd = pd.read_csv(RESULTS / "sampling_bounds.csv").set_index(["manual_station", "window"])
nm = {"Banaras Hindu University": "BHU", "Jawahar Nagar": "Jawahar Nagar", "Saket Nagar": "Saket Nagar", "Sigra": "Sigra"}
win = {k: k for k in s.window.unique()}
rows = []; below1 = 0
for r in s.itertuples():
    x = fd.loc[(r.manual_station, r.window)]
    assert abs(x.top9 - r.top9_days_each_month_bound) < 0.6 and abs(x.manual_mean - r.manual_mean) < 0.2, (r.manual_station, r.window)
    below1 += x.top1 < x.manual_mean
    rows.append([nm[r.manual_station], f"{r.continuous_station} ({r.distance_km} km)", win[r.window], f"{r.manual_mean:.0f}", f"{r.continuous_all_valid_days:.0f}",
                 f"{int(x.valid_days_total)}", f"{r.random9_p2_5:.0f} to {r.random9_p97_5:.0f}", f"{r.top9_days_each_month_bound:.0f}", f"{x.top1:.0f}"])
table("D2_sampling", ["Manual station", "Continuous station (distance)", "Period", "Manual mean", "Continuous, all valid days", "Valid days",
                      "Nine random days a month (95% of draws)", "Nine highest days a month", "Highest day a month"], rows)
print("D2 check: highest day a month below the manual mean in", int(below1), "of", len(rows))
bm = pd.read_csv(RESULTS / "sampling_bhu_months.csv")
mon = {"2023-11": "November 2023", "2023-12": "December 2023", "2024-01": "January 2024", "2024-02": "February 2024", "2024-03": "March 2024"}
rows = [[mon[r.month], f"{r.manual:.1f}", f"{r.continuous_mean:.1f}", int(r.valid_days), f"{r.top9:.0f}", f"{r.top4:.0f}", f"{r.highest_day:.0f}"] for r in bm.itertuples()]
rows.append(["Five months", f"{bm.manual.mean():.1f}", f"{bm.continuous_mean.mean():.1f}", int(bm.valid_days.sum()), f"{bm.top9.mean():.0f}", f"{bm.top4.mean():.0f}", f"{bm.highest_day.mean():.0f}"])
rows.append(["November to February", f"{bm.manual[:4].mean():.1f}", f"{bm.continuous_mean[:4].mean():.1f}", int(bm.valid_days[:4].sum()), f"{bm.top9[:4].mean():.0f}", f"{bm.top4[:4].mean():.0f}", f"{bm.highest_day[:4].mean():.0f}"])
assert (bm.manual > bm.highest_day).all()
table("D3_bhu_months", ["Month", "Manual", "Continuous, all days", "Valid days", "Nine highest days", "Four highest days", "Highest day"], rows)

# ------------------------------------------------------------------ E1 yearly series against the satellite estimate; E2 yearly difference; E3 matched chain
ns = pd.read_csv(RESULTS / "national_series.csv")
ns = ns[(ns.reference == "washu") & (ns.weights == "station-months")]
get = lambda g, sname: ns[(ns.group == g) & (ns.series == sname)].set_index("year")
rows = []
for y in range(2016, 2025):
    r = [y]
    for g in ["Uttar Pradesh", "Rest of India"]:
        a, b, c = get(g, "station").loc[y], get(g, "reference").loc[y], get(g, "ratio").loc[y]
        r += [int(c.stations_in_year), pc(a.pct), pc(b.pct), f"{pc(c.pct)} {iv(c.lo, c.hi)}"]
    rows.append(r)
table("E1_yearly", ["Year", "Uttar Pradesh: stations", "Station PM2.5", "Satellite PM2.5", "Station against satellite (95% interval)",
                    "Rest of India: stations", "Station PM2.5", "Satellite PM2.5", "Station against satellite (95% interval)"], rows)
gy = pd.read_csv(RESULTS / "gap_by_year.csv")
gy = gy[gy.stations == "All stations"]
rows = []
for y in [str(k) for k in range(2016, 2025)]:
    r = [y]
    for ref in ["washu", "cams", "ceres"]:
        x = gy[(gy.year == y) & (gy.reference == ref)].iloc[0]
        r += [f"{pc(x.up_pct)}, {pc(x.rest_pct)}", f"{m(x.diff_points, 0, plus=True)} {iv(x.lo, x.hi, '')}"]
    rows.append(r)
table("E2_gap", ["Year", "Satellite PM2.5: Uttar Pradesh, rest", "Difference, points", "CAMS PM2.5: Uttar Pradesh, rest", "Difference, points",
                 "Optical depth: Uttar Pradesh, rest", "Difference, points"], rows)
mt = pd.read_csv(RESULTS / "national_matched.csv")
mt = mt[mt.reference == "washu"]
rows = []
for y in range(2016, 2025):
    r = [y]
    for g in ["Uttar Pradesh", "Rest of India"]:
        x = mt[(mt.group == g) & (mt.year == y)].iloc[0]
        r += [int(x.stations_in_pair) if x.stations_in_pair else "", pc(x.chained_pct), pc(x.main_series_pct)]
    rows.append(r)
table("E3_matched", ["Year", "Uttar Pradesh: stations in both years", "Linked changes", "Main calculation",
                     "Rest of India: stations in both years", "Linked changes", "Main calculation"], rows)

# ------------------------------------------------------------------ F1 station table (166 eligible stations)
ss = pd.read_csv(RESULTS / "single_station.csv")
rt = pd.read_csv(RESULTS / "single_station_ratios.csv").set_index("station_id")
ss = ss.sort_values("washu_chg")
rows = [[short(r.station_name), r.state.replace("TamilNadu", "Tamil Nadu"), f"{rt.ratio_2016_19[r.station_id]:.2f}", f"{rt.ratio_2023_24[r.station_id]:.2f}",
         f"{m(r.washu_chg)}% ({m(r.washu_lo)} to {m(r.washu_hi)})", m(r.cams_chg) + "%", m(r.ceres_chg) + "%"] for r in ss.itertuples()]
table("F1_stations", ["Station", "State", "Ratio to satellite PM2.5, 2016–19", "Ratio, 2023–24", "Change against satellite PM2.5 (widened interval)",
                      "Against CAMS PM2.5", "Against optical depth"], rows)
a = ss[ss.station_id == "site_273"].iloc[0]
print("F1 check: n", len(ss), "UP", int((ss.state == "Uttar Pradesh").sum()), "Ardhali", round(a.washu_chg, 1), round(a.washu_lo), round(a.washu_hi),
      "wholly below", int((ss.washu_hi < 0).sum()), "above", int((ss.washu_lo > 0).sum()))

# ------------------------------------------------------------------ G record checks (tables lifted from the saved summaries)
def lift(path, header, name, head=None):
    hd, body = md_table(path, header)
    body = [[c.replace("-", "−") if re.fullmatch(r"[-+]?[\d.]+%?", c) else c for c in r] for r in body]
    table(name, head or hd, body)

FS = RESULTS / "ardhali_record_checks.md"
lift(FS, "## 2. Day-to-day coupling", "G1_coupling", head=["Station"] + [str(y) for y in range(2016, 2026)])
lift(FS, "## 3. Rain response", "G2_rain", head=["Station", "2016–19, wet against dry days, %", "Wet days", "2022–24, wet against dry days, %", "Wet days"])
lift(FS, "## 8c.", "G3_patna", head=["Pollutant", "Ratio to Patna, Jul 2020–Jun 2021", "Ratio to Patna, Aug 2021–Jul 2022", "Change"])
PS = RESULTS / "pm_ratio_within_month.md"
lift(PS, "## Ardhali Bazar, by month", "G4_ratio_months", head=["Year", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"])
lift(PS, "## All stations", "G5_ratio_national", head=["Year", "Stations", "Median station", "Stations with a fixed ratio", "Share"])

# ------------------------------------------------------------------ H1 input files, H2 result files
brk = lambda s: s.replace("_", "_​")   # zero-width space, so long names wrap in a narrow column
INPUTS = [("data/derived/stations.csv", "Continuous stations: name, city, state, coordinates and their source"),
          ("data/derived/continuous_station_month_pm25.csv", "Continuous PM2.5: station-month means and valid days"),
          ("data/derived/continuous_station_month_pm10.csv", "Continuous PM10: station-month means and valid days"),
          ("data/derived/satellite_pm25_station_month_v6.csv", "Satellite PM2.5 (V6.GL.03) at each station by month"),
          ("data/derived/satellite_pm25_station_month_v5.csv", "Satellite PM2.5 (V5.GL.06) at each station by month"),
          ("data/derived/manual_station_month.csv", "Manual stations: monthly values from the 84 reports"),
          ("data/derived/manual_city_centres.csv", "City-centre coordinates for cities with manual stations only"),
          ("data/derived/uppcb_reports.csv", "Addresses and checksums of the 84 reports"),
          ("data/raw/cams/eac4_monthly_india_2003-2025_pm.nc", "CAMS monthly PM2.5 and PM10 over India (downloaded)"),
          ("data/raw/weather/city_day.csv", "Daily weather by city (downloaded)"),
          ("data/raw/metar/VEBN.csv", "Varanasi airport reports (downloaded)")]
rows = []
for f, what in INPUTS:
    b = (ROOT / f).read_bytes()
    n = b.count(b"\n") - 1 if f.endswith(".csv") else ""
    rows.append([brk(Path(f).name), what, f"{n:,}" if n != "" else "", f"{len(b):,}", hashlib.sha256(b).hexdigest()[:12]])
table("H1_inputs", ["File", "Contents", "Rows", "Bytes", "SHA-256, first 12 characters"], rows)
RES = [("national_series.csv", "Yearly estimates by group and reference (Table 2, Table E1)", "11_national_series.py"),
       ("national_monthly.csv", "Monthly series (Figure 1)", "12_national_monthly.py"),
       ("national_steady.csv", "Steady-ratio test (Table 3)", "12_national_monthly.py"),
       ("national_matched.csv", "Linked year-on-year changes (Table E3)", "12_national_monthly.py"),
       ("single_station.csv", "Two-period change at 166 stations (Table F1, Figure 2)", "13_single_station.py"),
       ("gap_by_year.csv", "Yearly difference between groups (Table E2)", "14_gap_by_year.py"),
       ("pm25_base_periods.md", "PM2.5 by base period (Table A3)", "15_pm25_base_periods.py"),
       ("single_base_years.md", "Single base years (Table A3)", "16_single_base_years.py"),
       ("single_station_ratios.csv", "Ratios in each period (Table F1)", "17_single_station_checks.py"),
       ("two_period_variants.csv", "Other eligibility rules (Table A1)", "17_single_station_checks.py"),
       ("ardhali_yearly_ratio.csv", "Ardhali Bazar ratio by year (Figure 3)", "17_single_station_checks.py"),
       ("period_choice.csv", "Other base and end periods (Table A2)", "18_period_choice.py"),
       ("pm10_manual_continuous.csv", "Manual and continuous PM10 by year (Table 4)", "19_pm10_manual.py"),
       ("pm10_manual_checks.md", "Cities with both networks; flagged repeats kept (Table 4)", "20_pm10_manual_checks.py"),
       ("pm10_base_periods.md", "PM10 by base period (Table A3)", "21_pm10_base_periods.py"),
       ("pm25_pm10_ratio.md", "PM2.5 and PM10 at the same stations (Table A4)", "22_pm25_pm10_ratio.py"),
       ("step_tests/step_test_variants.csv", "Step test under eleven specifications (Table A6)", "23_step_tests.py, 24_step_summary.py"),
       ("step_tests/trend_vs_step.csv", "Trend and step models by station (Table A7)", "25_trend_and_dependence.py"),
       ("step_tests/residual_correlation_by_distance.csv", "Residual correlation by distance (Table A8)", "25_trend_and_dependence.py"),
       ("interval_coverage.csv", "Coverage and widening factor (Table A5)", "26_interval_coverage.py"),
       ("achieved_coverage.csv", "Coverage of the widened interval (Table A5)", "27_achieved_coverage.py"),
       ("weather_adjustment_results.csv", "Weather adjustment (Table B1, Figure 4)", "29_weather_adjustment.py"),
       ("airport_visibility_hourly.csv", "Airport visibility by hour (Table C1, Figure 4)", "30_airport_visibility.py"),
       ("airport_paired_pm25.csv", "Station PM2.5 in the airport hours (Table C1, Figure 4)", "31_airport_paired_pm25.py"),
       ("sampling_bounds.csv", "Sampling bounds (Table D2)", "32_sampling_bounds.py"),
       ("sampling_bhu_months.csv", "The university pair by month (Table D3)", "32_sampling_bounds.py"),
       ("ardhali_record_checks.md", "Checks on the Ardhali Bazar record (Tables G1 to G3)", "33_ardhali_record_checks.py"),
       ("pm_ratio_within_month.md", "Ratio of PM2.5 to PM10 within months (Tables G4, G5)", "34_pm_ratio_within_month.py")]
rows = []
for f, what, script in RES:
    for sc in script.split(", "):
        assert (ROOT / "code" / sc).exists(), sc
    rows.append([brk(Path(f).name), what, hashlib.sha256((RESULTS / f).read_bytes()).hexdigest()[:12], brk(script)])
table("H2_results", ["File in results/", "Contents", "SHA-256, first 12 characters", "Made by (code/)"], rows)
