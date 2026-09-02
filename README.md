# Remote sensing–based water volume and lithium potential of the Algero–Tunisian chotts

Code and derived data supporting:

> Hacini, M.E., Kallel, A., Kouadri, S., Hacini, M., Djidel, M., Khan, M.A., Alsubih, M., Islam, S.
> *Remote Sensing-Based Water Volume Estimation in Algero–Tunisian Chotts and Preliminary
> Assessment of Lithium Potential.* Submitted to *Pure and Applied Geophysics*.

**Scope.** This repository reproduces every number, table and figure derived from remote sensing in
the article. It does **not** contain a lithium resource or reserve estimate. The lithium figures are
a screening-level *potential*, `L = A × d × C`, computed from a remotely estimated water volume and
a brine concentration assumed from the literature. See *Known limitations* below before reusing
anything here.

Study areas: Chott Djerid and Chott Gharsa (Tunisia); Chott Melghigh and Chott Merouane (Algeria).
Period: January 2019 – December 2024, 550 cloud-screened scenes
(Djerid 144, Gharsa 136, Melghigh 135, Merouane 135).

## Data sources (all public)

| Purpose | Product | Earth Engine asset |
|---|---|---|
| Optical imagery (Gharsa, Melghigh, Merouane) | Landsat-8 OLI Collection 2 T1 L2 SR | `LANDSAT/LC08/C02/T1_L2` |
| Optical imagery (Djerid) | Sentinel-2 MSI L2A harmonised SR | `COPERNICUS/S2_SR_HARMONIZED` |
| Elevation | SRTM v3, 1 arc-second, void-filled (WGS84 / EGM96) | `USGS/SRTMGL1_003` |
| Rainfall | CHIRPS daily | `UCSB-CHG/CHIRPS/DAILY` |
| Air / dewpoint temperature | ERA5-Land monthly aggregated | `ECMWF/ERA5_LAND/MONTHLY_AGGR` |
| Potential evapotranspiration | TerraClimate | `IDAHO_EPSCOR/TERRACLIMATE` |
| Land-surface temperature | MODIS LST 8-day | `MODIS/061/MOD11A2` |

Analysis extents, per-chott sensor assignment and MNDWI thresholds:
`data/study_area_extents.geojson` (and `.csv`).

## Layout

```
scripts/
  calculateindices.py        GEE: image collection, cloud masking, spectral indices, MNDSI area export
  climate_extract_gee.py     GEE: CHIRPS / ERA5-Land / TerraClimate / MODIS extraction per chott
  climate_plot.py            Figure 9 (climate context) + enriched climate table
  water_volume_from_dem.py   DEM-based water depth and volume from a water mask  [see note below]
  analysis.py                MNDWI-based: Mann-Kendall + Theil-Sen trends, rainfall-lag and climate
                             correlations, lithium concentration sensitivity, depth statistics
  analysis_b.py              MNDSI-based: salt-affected area trends and water-salinity coupling
  make_figs.py               Figures 7 and 8 (requires analysis.py to have been run first)
data/
  Chott_*_water_volumes.csv          per-scene water area, volume, mean/max depth, surface elevation
  *_MNDSI.xlsx                       per-scene salt-affected area
  chotts_climate_2019_2024.csv       raw monthly climate extraction (GEE output)
  chotts_climate_enriched.csv        derived: temp_C, pet_mm, lst_C, rh_pct
  study_area_extents.geojson/.csv    analysis bounding boxes, sensor, MNDWI threshold per chott
  validation_confusion_matrix.csv    Table 3 (20 field points, coordinates withheld)
  trends.csv                         Table 4 (Mann-Kendall / Theil-Sen)
  correlations.csv                   Table 5 and Figure 7 (full lag spectrum)
  li_sensitivity.csv                 Table 8 and Figure 8
  salinity_trends_and_coupling.csv   MNDSI trends and water-salinity correlations
  depth_stats.csv                    per-chott depth and boundary-elevation statistics
  summary.json                       scene counts, date spans, thresholds, seasonal mean volumes
outputs/                             created by the scripts; git-ignored
```

## Reproducing the analysis

```bash
pip install -r requirements.txt
python scripts/analysis.py       # -> outputs/trends.csv, correlations.csv, li_sensitivity.csv,
                                 #    depth_stats.csv, summary.json
python scripts/analysis_b.py     # -> outputs/salinity_trends_and_coupling.csv
python scripts/make_figs.py      # -> outputs/Fig07..., Fig08...   (run analysis.py first)
python scripts/climate_plot.py   # -> outputs/Fig09..., chotts_climate_enriched.csv
```

Scripts resolve their paths relative to the repository root, so they can be run from anywhere.
The files written to `outputs/` are byte-identical to the corresponding files shipped in `data/`.

The two Earth Engine scripts (`calculateindices.py`, `climate_extract_gee.py`) require an
authenticated Earth Engine account (`earthengine authenticate`). `calculateindices.py` was written
for Google Colab and mounts Google Drive; adapt the two output paths to run it elsewhere.

## Method summary

1. **Water mask** — MNDWI = (Green − SWIR1)/(Green + SWIR1), thresholded at **−0.10** for Landsat-8
   and **−0.11** for Sentinel-2. Thresholds were calibrated separately per sensor because the SWIR1
   spectral response functions and spatial resolutions differ. Indices are computed at each sensor's
   native SWIR resolution (30 m Landsat-8, 20 m Sentinel-2).
2. **Salt-affected surface** — MNDSI = (SWIR1 − SWIR2)/(SWIR1 + SWIR2).
3. **Depth and volume** (`water_volume_from_dem.py`) — the mask is reprojected onto the SRTM grid and
   labelled into connected components; components smaller than 10 pixels are discarded; each
   component's boundary is `mask XOR binary_erosion(mask)`; the water surface of that component is
   assumed horizontal at the **mean elevation of its own boundary pixels**; depth is the positive
   difference with the DEM; volume is depth integrated over area. DEM nodata pixels are excluded.
4. **Lithium potential** — `L (t) = A (m²) × d (m) × C (mg L⁻¹) × 1e−6`, with C = 44.1 mg L⁻¹
   (Tunisian chotts) and 66 mg L⁻¹ (Algerian chotts) taken from the literature and tested over a
   0.5×–2.0× range.

## Not included

- **The water-mask export step.** `calculateindices.py` computes the spectral indices and exports
  MNDSI areas. The code that thresholded MNDWI and wrote the per-date water-mask rasters consumed by
  `water_volume_from_dem.py` is not included; the masks themselves are large rasters and are
  available from the corresponding author on request.
- **Field-validation coordinates**, which identify privately held land parcels. The aggregate
  confusion matrix is provided in `data/validation_confusion_matrix.csv`; coordinates are available
  from the corresponding author on reasonable request.

## Known limitations (please read before reusing these data)

- **No independent validation of depth or volume.** No bathymetric or GNSS-levelled depth data exist
  for these basins. Depth and water extent both derive from the same mask and DEM, so a systematic
  mask error propagates into the volume twice, in the same direction.
- **No sink filling was applied to the DEM.** Per-scene *maximum* depths of 26–37 m in
  `Chott_*_water_volumes.csv` are residual SRTM pits, not real depths. Use `Avg_Depth_m`; do not use
  `Max_Depth_m` as a physical quantity.
- **Lithium concentration is assumed, not measured.** Both literature values refer to near-surface
  brine from a small number of sampling points, with no published variance, and are applied
  uniformly across basins, seasons and years. Combined relative uncertainty on the lithium potential
  ranges from about ±66 % to more than ±300 %.
- **The MNDSI layer is unvalidated.** No electrical-conductivity or soil-salinity measurements are
  available. It maps surfaces whose SWIR response is consistent with salt enrichment; it quantifies
  no degree of salinity.
- **Six years is a short record.** No trend in annual mean water area or volume is significant at
  α = 0.05 for any chott. The significant results in `trends.csv` come from the monthly series, where
  the power derives from the number of observations rather than the length of the record.
- **The two sensors are not interchangeable.** Each basin's time series is single-sensor, so absolute
  areas are not strictly comparable between Chott Djerid (Sentinel-2) and the three Landsat-8 basins.

## Licence

MIT for the code, CC BY 4.0 for the derived data in `data/`. Source satellite products retain their
own terms. See `LICENSE`.

## Contact

Mohammed Elelmi Hacini — mohammed-elelmi.hacini@enis.tn
