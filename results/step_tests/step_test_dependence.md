# Step versus gradual trend, and dependence between stations

## Which model fits best (BIC; the break date counts as an extra parameter)

| group | stations | season only | trend | step | trend+step | median trend %/yr |
|---|---|---|---|---|---|---|
| Uttar Pradesh | 21 | 0 | 2 | 5 | 14 | -4.6 |
| Outside UP | 154 | 29 | 8 | 55 | 62 | -0.5 |
| UP, downward detections in the main specification | 9 | 0 | 1 | 3 | 5 | -7.8 |
| Outside UP, detections in the main specification | 43 | 0 | 3 | 30 | 10 | -3.5 |

Ardhali Bazar: best model trend+step; trend-only slope -14.5%/yr; BIC(step) - BIC(trend) = -11.9 (negative favours a step).
UP stations with a negative trend-only slope: 18 of 21; outside UP: 84 of 154.

## Correlation of step-model residuals between station pairs (monthly, >= 36 shared months)

| distance km | pairs | median r | mean r |
|---|---|---|---|
| 0-25 | 962 | 0.33 | 0.33 |
| 25-50 | 550 | 0.24 | 0.23 |
| 50-150 | 1085 | 0.17 | 0.17 |
| 150-300 | 1381 | 0.17 | 0.16 |
| 300-600 | 1914 | 0.12 | 0.12 |
| 600+ | 9198 | 0.03 | 0.03 |
| same city (any distance) | 744 | 0.36 | 0.36 |

## Global-null simulation with dependence kept (synchronized calendar blocks, 999 draws)

Each draw resamples the same calendar blocks for all 175 stations, so cross-station correlation of the residuals is kept.
Stations with p <= 0.05 per draw: mean 8.6 (independence: 8.8); SD 6.5 (independence: 2.9); 95th percentile 21; max 52.
Variance inflation of the count: 5.1; implied effective number of independent stations for counting purposes: 34 of 175.
BH 5% rejections per null draw: mean 0.05; share of draws with any rejection 0.006; 99th percentile 0; max 22.
Observed (sync p-values): 51 BH rejections; 72 with p <= 0.05.

## City clusters

UP: 11 cities with tested stations; 6 have at least one downward detection, 0 an upward one.
Outside UP: 85 cities; 23 with a downward detection, 13 with an upward one.
UP cities: Lucknow 3/4; Noida 0/4; Ghaziabad 2/3; Greater Noida 0/2; Meerut 0/2; Agra 1/1; Bulandshahr 1/1; Hapur 0/1; Kanpur 1/1; Muzaffarnagar 0/1; Varanasi 1/1

City-level permutation (cities as units, 96 cities, UP label shuffled across cities, 20,000 draws; exploratory, the UP contrast was chosen after seeing the data):
- net detection share (up minus down, per city): UP mean -0.49, others -0.08; one-sided p = 0.0119
- mean step estimate per city: UP -26.0%, others -2.2%; one-sided p = 0.0025
