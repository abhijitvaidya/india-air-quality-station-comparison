# Step test: detections under each specification (satellite PM2.5, V6.GL.03)

Counts: BH = Benjamini-Hochberg 5%; BY = Benjamini-Yekutieli 5% (valid under any dependence).

| variant | tested | bh_5pct | by_5pct | up_down | up_up | up_point_down | out_down | out_up | ardhali_step_pct | ardhali_month | ardhali_p | ardhali_bh | ardhali_rank_most_negative | description |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| orig | 175 | 52 | 27 | 9 | 0 | 18 | 27 | 16 | -58.3 | 202107 | 0.001 | True | 1 | Main specification: 4 season dummies, step null, 6-observation circular blocks over retained rows, trim 6 |
| trim24 | 175 | 45 | 21 | 11 | 0 | 19 | 20 | 14 | -58.3 | 202107 | 0.001 | True | 1 | Minimum segment 24 observations |
| calendar | 175 | 53 | 28 | 10 | 0 | 18 | 27 | 16 | -58.3 | 202107 | 0.001 | True | 1 | Calendar-month blocks of 6 (missing months kept missing) |
| sync | 175 | 51 | 27 | 9 | 0 | 18 | 27 | 15 | -58.3 | 202107 | 0.001 | True | 1 | Calendar blocks of 6 on common 2016-2024 grid, same block dates for all stations |
| month12 | 175 | 49 | 24 | 9 | 0 | 19 | 25 | 15 | -58.3 | 202107 | 0.001 | True | 1 | 12 month-of-year dummies instead of 4 seasons |
| trend | 175 | 64 | 0 | 2 | 5 | 9 | 34 | 23 | -43.8 | 202107 | 0.013 | True | 11 | Linear trend in null and alternative (step beyond trend) |
| drop_2020 | 175 | 60 | 30 | 10 | 0 | 19 | 33 | 17 | -61.0 | 202107 | 0.001 | True | 1 | Exclude all of 2020 |
| covid_dummy | 175 | 55 | 30 | 10 | 0 | 19 | 31 | 14 | -59.0 | 202107 | 0.001 | True | 1 | Mar-Jun 2020 dummy in both models |
| resid_null | 175 | 7 | 0 | 3 | 0 | 18 | 4 | 0 | -58.3 | 202107 | 0.001 | True | 1 | Bootstrap residuals from the no-step fit (not centred on the step model) |
| combined | 175 | 63 | 32 | 11 | 0 | 17 | 32 | 20 | -59.0 | 202107 | 0.001 | True | 1 | Calendar blocks 12, trim 12, exclude Mar-Jun 2020 |
| combined_trend | 175 | 60 | 0 | 2 | 5 | 10 | 28 | 25 | -46.2 | 202107 | 0.004 | True | 5 | As combined, plus linear trend (step beyond trend) |

## Overlap with the main specification

| variant | detected in both | only main specification | only variant |
|---|---|---|---|
| orig | 52 | 0 | 0 |
| trim24 | 39 | 13 | 6 |
| calendar | 52 | 0 | 1 |
| sync | 51 | 1 | 0 |
| month12 | 48 | 4 | 1 |
| trend | 30 | 22 | 34 |
| drop_2020 | 48 | 4 | 12 |
| covid_dummy | 50 | 2 | 5 |
| resid_null | 7 | 45 | 0 |
| combined | 47 | 5 | 16 |
| combined_trend | 27 | 25 | 33 |

## Detected under every specification except resid_null (10 specifications): 12 stations

| station | state | step %, main specification | min step % | max step % |
|---|---|---|---|---|
| Ardhali Bazar, Varanasi - UPPCB | Uttar Pradesh | -58 | -61 | -44 |
| Rabindra Bharati University, Kolkata - WBPCB | West Bengal | -45 | -46 | 77 |
| Khadakpada, Kalyan - MPCB | Maharashtra | -42 | -43 | -35 |
| Vikas Sadan, Gurugram - HSPCB | Haryana | -39 | -40 | -28 |
| DTU, Delhi - CPCB | Delhi | -33 | -34 | -33 |
| Sanjay Nagar, Ghaziabad - UPPCB | Uttar Pradesh | -30 | -31 | 52 |
| Sector-12, Karnal - HSPCB | Haryana | -30 | -33 | 35 |
| Borivali East, Mumbai - IITM | Maharashtra | -29 | -37 | 78 |
| Patti Mehar Ambala | Haryana | -28 | -36 | -25 |
| GVM Corporation, Visakhapatnam - APPCB | Andhra Pradesh | 20 | 18 | 27 |
| Asansol Court Area, Asansol - WBPCB | West Bengal | 44 | 38 | 69 |
| Samanpura, Patna - BSPCB | Bihar | 68 | 50 | 68 |

Robust set: UP 2 down / 0 up; outside UP 7 down / 3 up.

## Summary

Stations: 175 (21 in Uttar Pradesh).
Detected under all eight specifications without a trend: 31; Uttar Pradesh 8 down, 0 up; elsewhere 14 down, 9 up.
trend, Uttar Pradesh detections: Ardhali Bazar, Varanasi - UPPCB 202107 -44%; Lalbagh, Lucknow - CPCB 202103 -31%; Sector - 62, Noida - IMD 201806 +39%; Sanjay Nagar, Ghaziabad - UPPCB 202312 +52%; Sector-1, Noida - UPPCB 202402 +54%; Knowledge Park - V, Greater Noida - UPPCB 202308 +58%; Indirapuram, Ghaziabad - UPPCB 202404 +95%
combined_trend, Uttar Pradesh detections: Ardhali Bazar, Varanasi - UPPCB 202107 -46%; Loni, Ghaziabad - UPPCB 202209 -30%; Knowledge Park - III, Greater Noida - UPPCB 202211 +22%; Sanjay Nagar, Ghaziabad - UPPCB 202312 +52%; Sector-1, Noida - UPPCB 202401 +52%; Knowledge Park - V, Greater Noida - UPPCB 202308 +63%; Indirapuram, Ghaziabad - UPPCB 202401 +76%
resid_null, Uttar Pradesh detections: Ardhali Bazar, Varanasi - UPPCB 202107 -58%; Sanjay Palace, Agra - UPPCB 202207 -51%; Talkatora District Industries Center, Lucknow - CPCB 202105 -29%
resid_null, all detections: DTU, Delhi - CPCB 202209 -33%; Ardhali Bazar, Varanasi - UPPCB 202107 -58%; Rabindra Bharati University, Kolkata - WBPCB 202301 -45%; Talkatora District Industries Center, Lucknow - CPCB 202105 -29%; Manali, Chennai - CPCB 202012 -43%; Sanjay Palace, Agra - UPPCB 202207 -51%; Maninagar, Ahmedabad - GPCB 201911 -32%

## Earlier satellite version, V5.GL.06

orig: 164 stations; Uttar Pradesh 8 down, 0 up; Ardhali Bazar -52% in 202107.
trend: 164 stations; Uttar Pradesh 1 down, 5 up; Ardhali Bazar -41% in 202107.
Correlation of step estimates, V6 and V5, main specification: 0.90 over 164 stations.
