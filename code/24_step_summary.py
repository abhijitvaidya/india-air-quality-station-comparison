#!/usr/bin/env python3
"""Counts of step detections under each specification (Table A6), stations detected under every specification, and the
comparison with the earlier satellite version V5.GL.06 (Appendix A8).

Run after 23_step_tests.py:  python code/24_step_summary.py
Outputs: results/step_tests/step_test_variants.csv/.md, results/step_tests/step_detection_matrix.csv
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import *
import os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import paired, bh, by, SEASONS, HERE, DS
R = str(STEP)
DESC = {
 'orig': 'Main specification: 4 season dummies, step null, 6-observation circular blocks over retained rows, trim 6',
 'L3': 'Blocks of 3 retained observations', 'L12': 'Blocks of 12 retained observations', 'L24': 'Blocks of 24 retained observations',
 'trim12': 'Minimum segment 12 observations', 'trim24': 'Minimum segment 24 observations',
 'calendar': 'Calendar-month blocks of 6 (missing months kept missing)', 'calendar_L12': 'Calendar-month blocks of 12',
 'sync': 'Calendar blocks of 6 on common 2016-2024 grid, same block dates for all stations',
 'month12': '12 month-of-year dummies instead of 4 seasons', 'month12_calendar_L6': '12 month dummies, calendar blocks of 6',
 'month12_calendar_L12': '12 month dummies, calendar blocks of 12',
 'trend': 'Linear trend in null and alternative (step beyond trend)',
 'trend_calendar_L12_trim12': 'Trend + calendar blocks 12 + trim 12',
 'drop_lockdown': 'Exclude Mar-Jun 2020', 'drop_2020': 'Exclude all of 2020', 'covid_dummy': 'Mar-Jun 2020 dummy in both models',
 'resid_null': 'Bootstrap residuals from the no-step fit (not centred on the step model)',
 'combined': 'Calendar blocks 12, trim 12, exclude Mar-Jun 2020',
 'combined_month12': 'As combined, with 12 month dummies',
 'combined_trend': 'As combined, plus linear trend (step beyond trend)',
}
rows = []
for v in DESC:
    f = f'{R}/break_{v}_v6.csv'
    if not os.path.exists(f): continue
    r = pd.read_csv(f, dtype={'break_month': str})
    up = r.state.eq('Uttar Pradesh')
    a = r[r.station_id == 'site_273'].iloc[0]
    rows.append({'variant': v, 'description': DESC[v], 'tested': len(r), 'bh_5pct': int(r.bh_5pct.sum()), 'by_5pct': int(r.by_5pct.sum()),
                 'p_le_0.05': int((r.p <= 0.05).sum()),
                 'up_down': int((up & r.bh_5pct & (r.shift_pct < 0)).sum()), 'up_up': int((up & r.bh_5pct & (r.shift_pct > 0)).sum()),
                 'up_n': int(up.sum()), 'up_point_down': int((up & (r.shift_pct < 0)).sum()),
                 'out_down': int((~up & r.bh_5pct & (r.shift_pct < 0)).sum()), 'out_up': int((~up & r.bh_5pct & (r.shift_pct > 0)).sum()),
                 'out_n': int((~up).sum()),
                 'up_down_cities': r[up & r.bh_5pct & (r.shift_pct < 0)].city.nunique(),
                 'out_down_cities': r[~up & r.bh_5pct & (r.shift_pct < 0)].city.nunique(),
                 'out_up_cities': r[~up & r.bh_5pct & (r.shift_pct > 0)].city.nunique(),
                 'ardhali_step_pct': round(a.shift_pct, 1), 'ardhali_month': a.break_month, 'ardhali_p': a.p, 'ardhali_bh': bool(a.bh_5pct),
                 'ardhali_rank_most_negative': int((r.shift_pct < a.shift_pct).sum()) + 1})
t = pd.DataFrame(rows)
t.to_csv(f'{R}/step_test_variants.csv', index=False)
cols = ['variant', 'tested', 'bh_5pct', 'by_5pct', 'up_down', 'up_up', 'up_point_down', 'out_down', 'out_up', 'ardhali_step_pct', 'ardhali_month', 'ardhali_p', 'ardhali_bh', 'ardhali_rank_most_negative']
md = ['# Step test: detections under each specification (satellite PM2.5, V6.GL.03)', '',
      'Counts: BH = Benjamini-Hochberg 5%; BY = Benjamini-Yekutieli 5% (valid under any dependence).', '',
      '| ' + ' | '.join(cols) + ' | description |', '|' + '---|' * (len(cols) + 1)]
for _, x in t.iterrows():
    md.append('| ' + ' | '.join(str(x[c]) for c in cols) + f' | {x.description} |')

# agreement of detections with the main specification
o = pd.read_csv(f'{R}/break_orig_v6.csv').set_index('station_id')
md += ['', '## Overlap with the main specification', '', '| variant | detected in both | only main specification | only variant |', '|---|---|---|---|']
for v in t.variant:
    r = pd.read_csv(f'{R}/break_{v}_v6.csv').set_index('station_id')
    a, b = o.bh_5pct, r.bh_5pct.reindex(o.index).fillna(False)
    md.append(f'| {v} | {int((a & b).sum())} | {int((a & ~b).sum())} | {int((~a & b).sum())} |')
# stations robust across all main variants (exclude resid_null)
keep = [v for v in t.variant if v != 'resid_null']
det = pd.concat([pd.read_csv(f'{R}/break_{v}_v6.csv').set_index('station_id').bh_5pct.rename(v) for v in keep], axis=1).fillna(False)
sh = pd.concat([pd.read_csv(f'{R}/break_{v}_v6.csv').set_index('station_id').shift_pct.rename(v) for v in keep], axis=1)
robust = det.all(1)
info = o[['station_name', 'state', 'shift_pct']]
rb = info[robust].sort_values('shift_pct')
md += ['', f'## Detected under every specification except resid_null ({len(keep)} specifications): {int(robust.sum())} stations', '',
       '| station | state | step %, main specification | min step % | max step % |', '|---|---|---|---|---|']
for sid, x in rb.iterrows():
    md.append(f'| {x.station_name} | {x.state} | {x.shift_pct:.0f} | {sh.loc[sid].min():.0f} | {sh.loc[sid].max():.0f} |')
up = o.state.eq('Uttar Pradesh')
md += ['', f'Robust set: UP {int((robust & up & (o.shift_pct < 0)).sum())} down / {int((robust & up & (o.shift_pct > 0)).sum())} up; '
       f'outside UP {int((robust & ~up & (o.shift_pct < 0)).sum())} down / {int((robust & ~up & (o.shift_pct > 0)).sum())} up.']
det.join(info).to_csv(f'{R}/step_detection_matrix.csv')
open(f'{R}/step_test_variants.md', 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))

# Detections under every specification without a trend that uses residuals of the step model, and detections in
# Uttar Pradesh with a trend allowed (Appendix A6)
rd = lambda v, ver='v6': pd.read_csv(f'{R}/break_{v}_{ver}.csv', dtype={'break_month': str}).set_index('station_id')
NT = ['orig', 'month12', 'calendar', 'sync', 'trim24', 'drop_2020', 'covid_dummy', 'combined']
o0 = rd('orig'); all8 = pd.concat([rd(v).bh_5pct.rename(v) for v in NT], axis=1).all(1); u = o0.state == 'Uttar Pradesh'
X = ['', '## Summary', '', f'Stations: {len(o0)} ({int(u.sum())} in Uttar Pradesh).',
     f'Detected under all eight specifications without a trend: {int(all8.sum())}; Uttar Pradesh {int((all8 & u & (o0.shift_pct < 0)).sum())} down, '
     f'{int((all8 & u & (o0.shift_pct > 0)).sum())} up; elsewhere {int((all8 & ~u & (o0.shift_pct < 0)).sum())} down, {int((all8 & ~u & (o0.shift_pct > 0)).sum())} up.']
for v in ['trend', 'combined_trend', 'resid_null']:
    t = rd(v); d = t[t.bh_5pct & (t.state == 'Uttar Pradesh')].sort_values('shift_pct')
    X.append(f'{v}, Uttar Pradesh detections: ' + '; '.join(f'{r.station_name} {r.break_month} {r.shift_pct:+.0f}%' for r in d.itertuples()))
t = rd('resid_null'); X.append('resid_null, all detections: ' + '; '.join(f'{r.station_name} {r.break_month} {r.shift_pct:+.0f}%' for r in t[t.bh_5pct].itertuples()))
# earlier satellite version
if os.path.exists(f'{R}/break_orig_v5.csv'):
    X += ['', '## Earlier satellite version, V5.GL.06', '']
    for v in ['orig', 'trend']:
        t = rd(v, 'v5'); uu = t.state == 'Uttar Pradesh'; a = t.loc['site_273']
        X.append(f'{v}: {len(t)} stations; Uttar Pradesh {int((uu & t.bh_5pct & (t.shift_pct < 0)).sum())} down, {int((uu & t.bh_5pct & (t.shift_pct > 0)).sum())} up; '
                 f'Ardhali Bazar {a.shift_pct:+.0f}% in {a.break_month}.')
    b = rd('orig').join(rd('orig', 'v5'), rsuffix='_v5', how='inner')
    X.append(f'Correlation of step estimates, V6 and V5, main specification: {b.shift_pct.corr(b.shift_pct_v5):.2f} over {len(b)} stations.')
open(f'{R}/step_test_variants.md', 'a').write('\n'.join(X) + '\n'); print('\n'.join(X))
