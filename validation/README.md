# SIAC validation cases

SIAC output compared with RadCalNet, AERONET and USGS Collection 2 surface
reflectance. Each case folder holds the numbers (`*.csv`) and the exact run
identity (`case.json`: product IDs, SIAC commit and config, region, reference
record files and times).

```python
import pandas as pd
df = pd.read_csv("validation/cases/landsat9_rvus_2025-06-16/results.csv", index_col="band")
```

## Cases

Differences are product minus reference, in reflectance (or AOD) units.

| Case | Compared | Headline | Caveat |
|---|---|---|---|
| [landsat8_rvus_2025-07-10](cases/landsat8_rvus_2025-07-10) | L8 OLI: SIAC v2, SIAC v1, USGS SR vs RadCalNet RVUS (same day) | All three 0.007-0.025 below RadCalNet in bands 1-4, within 0.006 in bands 5-6, +0.021 to +0.027 in band 7 | One scene |
| [landsat8_rvus_2024-11-28](cases/landsat8_rvus_2024-11-28) | L8 OLI: same three vs RadCalNet RVUS records 4 days before and 2 days after | All three 0.009-0.049 below the closer (30 Nov) record in every band | No same-day reference |
| [landsat9_rvus_2025-06-16](cases/landsat9_rvus_2025-06-16) | L9 OLI-2: SIAC v2, USGS SR vs RadCalNet surface; Level-1 TOA vs RadCalNet TOA | Level-1 TOA is already 0.018-0.023 low in bands 1-4 and +0.019 in band 7; SIAC v2 surface follows that pattern. USGS SR is -0.106 / -0.067 in bands 1 / 2, not present at TOA | RadCalNet TOA is modelled, not measured |
| [siac_v1_sentinel2_radcalnet](cases/siac_v1_sentinel2_radcalnet) | S2A/B/C: SIAC v1 BOA and its stated uncertainty vs RadCalNet, 5 sites, 113 scenes | 1051 band samples: median abs. residual 0.009, 78 % within SIAC's 1 sigma (68 % expected); per site and band in `coverage_by_site_band.csv` | BTCN and the BSCN 20/60 m bands left out: window larger than the target |
| [siac_v1_sentinel2_aeronet_aod](cases/siac_v1_sentinel2_aeronet_aod) | S2A/B/C: SIAC v1 AOD550 vs AERONET, 7 stations, 87 matchups | RMSE 0.127, bias -0.036 overall; RMSE 0.061 below AOD 0.2; retrieval within RMSE 0.020 of the CAMS prior it started from | Stations weighted to heavy aerosol |

A difference against RadCalNet or AERONET is not by itself a SIAC error: the
reference, sensor calibration, spectral response, spatial support and timing all
contribute. Single scenes are not accuracy figures.

## How the numbers were made

- **Landsat region.** 1 km x 1 km north-up square on the RadCalNet site
  coordinates; each 30 m pixel weighted by its fractional overlap.
- **Sentinel-2 region.** Median over a 21 x 21 native-pixel window (0.21 / 0.42 /
  1.26 km for 10 / 20 / 60 m bands).
- **RadCalNet to band values.** 10 nm spectrum weighted by the band's published
  relative spectral response (USGS for Landsat OLI/OLI-2, ESA v4.0 for
  Sentinel-2 MSI). Per-wavelength uncertainty treated as fully correlated within
  a band (upper bound) and as k=1, since RadCalNet publishes no coverage factor.
- `*_roi_sd` is spatial variation inside the region, not an uncertainty.
- The Landsat cases use the Level-1 preprocessor on branch
  `feature/landsat-c2-l1` of this fork (commit 6f4028e; the runs used
  identical code as 70b3570, before `feat/refactor` was rewritten); SIAC v2 config in
  `cases/config/siac_v2_landsat.toml`.

`tools/fetch_radcalnet.py` re-downloads the RadCalNet records a case used and
recomputes its reference columns (uses spectrAccess).
