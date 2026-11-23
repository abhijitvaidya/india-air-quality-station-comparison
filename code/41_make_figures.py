#!/usr/bin/env python3
"""Figures 1 to 4 of the paper, from the result files.

Run after the analysis scripts:  python code/41_make_figures.py
Writes figures/fig1..fig4 as PNG (300 dpi) and PDF.

Colour by what a series is, the same in every figure:
  Uttar Pradesh and its stations orange, rest of India blue, references and other stations grey,
  the three later Varanasi stations magenta. Each series also has its own marker or line style,
  so the figures read in greyscale.
"""
import os, re, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
HERE = FIGURES
M = RESULTS
DS = DERIVED

BLUE, ORANGE, MAGENTA = "#2a78d6", "#eb6834", "#e87ba4"
INK, INK2, GRID, GREY = "#0b0b0b", "#52514e", "#e6e5e1", "#a9a8a2"
plt.rcParams.update({
    "font.size": 9, "font.family": "DejaVu Sans", "text.color": INK, "axes.labelcolor": INK2,
    "axes.edgecolor": GRID, "axes.linewidth": 0.8, "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True, "legend.frameon": False,
    "legend.fontsize": 8, "figure.dpi": 100, "savefig.dpi": 300, "pdf.fonttype": 42, "axes.titlesize": 9,
    "axes.titleweight": "bold", "axes.titlelocation": "left"})
PF = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:+.0f}%".replace("-", "−") if v else "0")


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


# ------------------------------------------------------------------ Figure 1: monthly series by group
def fig1():
    s = pd.read_csv(M / "national_monthly.csv", dtype={"ym": str})
    s["t"] = s.ym.str[:4].astype(int) + (s.ym.str[4:].astype(int) - 0.5) / 12
    fig = plt.figure(figsize=(7.4, 3.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[4.2, 1], hspace=0.1, wspace=0.1)
    ax = [fig.add_subplot(gs[0, i]) for i in range(3)]
    bx = [fig.add_subplot(gs[1, i]) for i in range(3)]

    def get(g, what, ref="washu"):
        return s[(s.group == g) & (s.reference == ref) & (s.series == what)].sort_values("t")

    def line(a, g, what, col, ls="-", label=None, band=False, dots=True, lw=1.9):
        x = get(g, what)
        if dots: a.plot(x.t, x.pct, ls="", marker="o", ms=1.6, color=col, alpha=0.35, mew=0)
        if band: a.fill_between(x.t, x.lo_12m, x.hi_12m, color=col, alpha=0.14, lw=0)
        a.plot(x.t, x.pct_12m, color=col, ls=ls, lw=lw, label=label)

    for a, b, g, col, t in [(ax[0], bx[0], "Uttar Pradesh", ORANGE, "a  Uttar Pradesh"), (ax[1], bx[1], "Rest of India", BLUE, "b  Rest of India")]:
        line(a, g, "station", col, label="Station PM2.5")
        line(a, g, "reference", INK2, ls=(0, (4, 2)), label="Satellite PM2.5", dots=False, lw=1.4)
        a.set_title(t); a.legend(loc="upper right", fontsize=7)
        x = get(g, "ratio")
        b.fill_between(x.t, 0, x.stations, step="mid", color=col, alpha=0.35, lw=0); b.plot(x.t, x.stations, color=col, lw=1, drawstyle="steps-mid")
        b.set_ylim(0, x.stations.max() * 1.45); b.set_yticks([])
        b.annotate(f"{int(x.stations.iloc[0])} stations", (x.t.iloc[0], x.stations.iloc[0]), xytext=(2, 5), textcoords="offset points", fontsize=7, color=INK2)
        b.annotate(f"{int(x.stations.iloc[-1])}", (x.t.iloc[-1], x.stations.iloc[-1]), xytext=(-2, 3), textcoords="offset points", fontsize=7, color=INK2, ha="right")
    line(ax[2], "Uttar Pradesh", "ratio", ORANGE, label="Uttar Pradesh", band=True)
    line(ax[2], "Rest of India", "ratio", BLUE, label="Rest of India", band=True)
    line(ax[2], "Plain outside Uttar Pradesh", "ratio", INK2, ls=":", label="Plain outside Uttar Pradesh", dots=False, lw=1.5)
    ax[2].set_title("c  Station / satellite ratio"); ax[2].legend(loc="upper right", fontsize=7)
    for i, a in enumerate(ax):
        a.axhline(0, color=GREY, lw=0.8); a.set_xlim(2016, 2025); a.set_ylim(-70, 62); a.set_xticks(range(2017, 2025, 2)); a.tick_params(labelbottom=(i == 2))
        a.yaxis.set_major_formatter(PF)
        if i: a.tick_params(labelleft=False)
    for b in bx[:2]:
        b.set_xlim(2016, 2025); b.set_xticks(range(2017, 2025, 2)); b.grid(False)
    bx[2].axis("off")
    bx[2].text(0.0, 0.3, "Dots: single months. Lines: 12-month mean.\nBands: 95% intervals, cities resampled.\nLower strips: stations with data.", fontsize=7, color=INK2, va="center", transform=bx[2].transAxes)
    ax[0].set_ylabel("Change from the same months\nof 2016–2019")
    save(fig, "fig1_monthly_groups")
    return s


# ------------------------------------------------------------------ Figure 2: 166 eligible stations, two references
def fig2():
    d = pd.read_csv(M / "single_station.csv")
    up = d.state == "Uttar Pradesh"
    fig, ax = plt.subplots(figsize=(6.3, 5.0))
    lim = (-90, 90)
    ax.plot(lim, lim, color=GREY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.axhline(0, color=GREY, lw=0.8, zorder=1); ax.axvline(0, color=GREY, lw=0.8, zorder=1)
    ax.scatter(d.washu_chg[~up], d.cams_chg[~up], s=20, color=GREY, edgecolor="white", lw=0.5, label=f"Stations elsewhere ({int((~up).sum())})", zorder=2)
    ax.scatter(d.washu_chg[up], d.cams_chg[up], s=30, color=ORANGE, edgecolor="white", lw=0.6, marker="D", label=f"Stations in Uttar Pradesh ({int(up.sum())})", zorder=3)
    lab = {"site_273": ("Ardhali Bazar, Varanasi", (14, -6), "left"), "site_307": ("Sanjay Palace, Agra", (-4, 34), "center")}
    for sid, (txt, off, ha) in lab.items():
        r = d[d.station_id == sid].iloc[0]
        ax.annotate(txt, (r.washu_chg, r.cams_chg), xytext=off, textcoords="offset points", ha=ha, va="center",
                    fontsize=8, color=INK, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6, shrinkA=0, shrinkB=3))
    ax.text(18, 21, "same change against both", fontsize=7.5, color=INK2, rotation=38, rotation_mode="anchor", ha="center", va="bottom", clip_on=True)
    XL, YL = (-75, 65), (-85, 65)
    ax.set_xlim(*XL); ax.set_ylim(*YL)
    ax.set_xlabel("Change in station / satellite PM2.5 ratio, 2016–2019 to 2023–2024")
    ax.set_ylabel("Change in station / CAMS PM2.5 ratio")
    ax.xaxis.set_major_formatter(PF); ax.yaxis.set_major_formatter(PF); ax.legend(loc="upper left", handletextpad=0.3)
    out = d[(d.washu_chg > XL[1]) | (d.cams_chg > YL[1]) | (d.washu_chg < XL[0]) | (d.cams_chg < YL[0])]
    for k, r in enumerate(out.itertuples()):   # beyond the axes: an arrow at the edge, with the values
        ax.annotate(f"{re.sub(r' - [^-]*$', '', r.station_name)}\n({r.washu_chg:+.0f}%, {r.cams_chg:+.0f}%)", xy=(64, 64), xytext=(44, 57 - 12 * k),
                    ha="center", va="center", fontsize=7.5, color=INK2, arrowprops=dict(arrowstyle="->", color=INK2, lw=0.8))
    save(fig, "fig2_stations_two_references")
    return d, out


# ------------------------------------------------------------------ Figure 3: Ardhali Bazar ratio over time
def fig3():
    g = pd.read_csv(DS / "continuous_station_month_pm25.csv")
    s = pd.read_csv(DS / "satellite_pm25_station_month_v6.csv")
    m = g[g.ndays >= 15].merge(s, on=["station_id", "ym"])
    m = m[(m.ym >= 201601) & (m.ym <= 202412)].copy()
    m["ratio"] = m["mean"] / m.sat_pm25
    m["t"] = pd.to_datetime(m.ym.astype(str), format="%Y%m") + pd.Timedelta(days=14)
    ay = pd.read_csv(RESULTS / "ardhali_yearly_ratio.csv")
    ann = {int(y): round(float(v), 2) for y, v in zip(ay.year, ay.ratio)}
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    new = {"site_5469": "Bhelupur", "site_5468": "IESD BHU", "site_5475": "Maldahiya"}
    for i, sid in enumerate(new):
        x = m[m.station_id == sid].sort_values("t")
        ax.plot(x.t, x.ratio, color=MAGENTA, lw=1.0, alpha=0.9, label="Three stations that began in June–July 2021" if i == 0 else None, zorder=2)
    a = m[m.station_id == "site_273"].sort_values("t")
    ax.plot(a.t, a.ratio, color=ORANGE, lw=1.0, alpha=0.55, zorder=3)
    ax.scatter(a.t, a.ratio, s=7, color=ORANGE, zorder=3, label="Ardhali Bazar, monthly")
    yrs = sorted(ann)
    ax.plot([pd.Timestamp(y, 7, 1) for y in yrs], [ann[y] for y in yrs], color=INK, lw=1.6, marker="o", ms=5,
            mfc="white", mew=1.4, label="Ardhali Bazar, yearly (seasons weighted equally)", zorder=4)
    for y in yrs:
        ax.annotate(f"{ann[y]:.2f}", (pd.Timestamp(y, 7, 1), ann[y]), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=7.5, color=INK)
    ax.axvline(pd.Timestamp(2021, 7, 1), color=INK2, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.text(pd.Timestamp(2021, 7, 20), 3.3, "New stations begin\n22 June and 8 July 2021", fontsize=7.5, color=INK2, va="top")
    ax.axhline(1, color=GREY, lw=0.8, zorder=1)
    ax.set_yscale("log"); ax.set_ylim(0.09, 3.6)
    ax.set_yticks([0.125, 0.25, 0.5, 1, 2]); ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2"]); ax.minorticks_off()
    ax.set_ylabel("Station PM2.5 / satellite PM2.5")
    ax.set_xlim(pd.Timestamp(2016, 1, 1), pd.Timestamp(2025, 1, 1))
    ax.legend(loc="lower left", ncol=1)
    save(fig, "fig3_ardhali_ratio")
    return ann


# ------------------------------------------------------------------ Figure 4: checks at Ardhali Bazar
def fig4():
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.0), gridspec_kw={"width_ratios": [1.15, 1.05, 1.05]})
    # (a) weather-adjusted change 2019 -> 2024
    w = pd.read_csv(RESULTS / "weather_adjustment_results.csv")
    w = w[w.spec.str.startswith("C ")].set_index("city").loc[["Varanasi", "Kanpur", "Lucknow", "Patna", "Delhi"]]
    ax = axes[0]; y = np.arange(len(w))[::-1]
    names = ["Varanasi\n(Ardhali Bazar)", "Kanpur", "Lucknow", "Patna", "Delhi"]
    for yi, (city, r) in zip(y, w.iterrows()):
        col = ORANGE if city == "Varanasi" else GREY
        ax.plot([r.weather_2024_lo, r.weather_2024_hi], [yi, yi], color=col, lw=2, solid_capstyle="round", zorder=2)
        ax.scatter([r.weather_2024], [yi], s=34, color=col, edgecolor="white", lw=0.8, zorder=3, marker="D" if city == "Varanasi" else "o")
        ax.annotate(f"{r.weather_2024:.0f}%".replace("-", "−"), (r.weather_2024_hi, yi), xytext=(5, 0), textcoords="offset points", va="center", fontsize=7.5, color=INK)
    ax.axvline(0, color=GREY, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7.5); ax.set_xlim(-85, 24); ax.set_xticks([-75, -50, -25, 0]); ax.xaxis.set_major_formatter(PF)
    ax.grid(axis="y", visible=False)
    ax.set_title("a  Weather-adjusted PM2.5,\n    2019 to 2024")
    # (b) airport visibility: share of selected daytime hours below 4 km (one observation per airport-hour)
    q = pd.read_csv(RESULTS / "airport_visibility_hourly.csv")
    ax = axes[1]
    for icao, name, col, mk, yy in [("VEPT", "Patna airport", GREY, "o", 21), ("VEBN", "Varanasi airport", INK, "s", 53)]:
        g = q[(q.icao == icao) & (q.year >= 2016) & (q.year <= 2024)].sort_values("year")
        ax.plot(g.year, 100 * g.hour_median_lt4000, color=col, lw=1.8, marker=mk, ms=4.5, mec="white", mew=0.6, zorder=3)
        ax.annotate(name, (2016.1, yy), ha="left", fontsize=7.5, color=INK if icao == "VEBN" else INK2)
    ax.set_ylim(0, 60); ax.set_xticks([2016, 2018, 2020, 2022, 2024])
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    ax.set_title("b  Dry daytime hours with\n    visibility below 4 km")
    # (c) Ardhali Bazar PM2.5 in the hours paired with the Varanasi airport reports
    p = pd.read_csv(RESULTS / "airport_paired_pm25.csv").sort_values("y")
    p = p[(p.y >= 2016) & (p.y <= 2024)]
    ax = axes[2]
    ax.plot(p.y, p.cpcb_pm_med, color=ORANGE, lw=1.8, marker="D", ms=4.5, mec="white", mew=0.6, zorder=3)
    for yy in (2019, 2024):
        v = float(p[p.y == yy].cpcb_pm_med.iloc[0])
        ax.annotate(f"{v:.0f}", (yy, v), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7.5, color=INK)
    ax.set_ylim(0, max(80, p.cpcb_pm_med.max() * 1.2)); ax.set_xticks([2016, 2018, 2020, 2022, 2024])
    ax.set_title("c  Ardhali Bazar PM2.5 in the\n    same hours, median, µg/m³")
    fig.tight_layout(w_pad=1.6)
    save(fig, "fig4_ardhali_checks")
    return w, q, p


if __name__ == "__main__":
    for old in list(HERE.glob("fig*.png")) + list(HERE.glob("fig*.pdf")):
        old.unlink()
    s = fig1(); d, out = fig2(); ann = fig3(); w, q, p = fig4()
    print("fig1 groups", sorted(s.group.unique()))
    print("fig2 n", len(d), "UP", int((d.state == "Uttar Pradesh").sum()), "x", d.washu_chg.min(), d.washu_chg.max(), "y", d.cams_chg.min(), d.cams_chg.max(), "beyond axes", len(out))
    print("fig3 annual", ann)
    print("fig4 paired PM2.5", dict(zip(p.y, p.cpcb_pm_med)), "| below 4 km VEBN 2019/2024", [round(100 * float(q[(q.icao == 'VEBN') & (q.year == yy)].hour_median_lt4000.iloc[0])) for yy in (2019, 2024)])
