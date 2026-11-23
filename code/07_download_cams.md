# CAMS EAC4 monthly fields

The CAMS reanalysis file is downloaded by hand from the Copernicus Atmosphere Data Store (free account needed):

1. Open the dataset "CAMS global reanalysis (EAC4) monthly averaged fields":
   https://ads.atmosphere.copernicus.eu/datasets/cams-global-reanalysis-eac4-monthly
2. Choose product type "Monthly mean", variables "Particulate matter d < 2.5 µm (PM2.5)" and
   "Particulate matter d < 10 µm (PM10)", years 2003 to 2025, all months, and the sub-region
   North 37.5, West 68.25, South 6, East 97.5. Download as GRIB.
3. Convert to netCDF (for example `cfgrib` with `xarray`: `xr.open_dataset('data.grib', engine='cfgrib').to_netcdf(...)`)
   and save as `data/raw/cams/eac4_monthly_india_2003-2025_pm.nc`.

The file used in the paper was retrieved on 17 November 2026 and has 276 months on a 43 × 40 grid (0.75 degree).
Its SHA-256 is listed in `data/README.md`.
