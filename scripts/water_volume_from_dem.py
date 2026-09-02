# -*- coding: utf-8 -*-
"""
Water volume estimation from an MNDWI water mask and the SRTM DEM.

This is the routine that produced data/Chott_*_water_volumes.csv, i.e. every
water volume in Table 6 of the article and, through them, every lithium figure
in Table 7. It is reproduced here in full because the DEM-based volume estimate
is the dominant source of uncertainty in the study and must be auditable.

METHOD (Section 3.3 of the article)
  1. The binary water mask for one acquisition date is reprojected onto the DEM
     grid and labelled into connected components (scipy.ndimage.label, default
     4-connectivity).
  2. Components smaller than MIN_PIXELS are discarded as classification noise.
  3. For each remaining component the boundary is obtained as the set difference
     between the component and its binary erosion (default 4-connectivity
     structuring element), i.e. the inner ring of pixels that touch a non-water
     pixel.
  4. The water surface of that component is assumed horizontal at the MEAN
     elevation of its boundary pixels.
  5. Depth is the positive difference between that surface elevation and the DEM;
     negative values are set to zero.
  6. Volume is the sum of depth x pixel area over the component.
  7. Per-date totals are the sum over all components.

IMPORTANT CAVEATS (stated in the article and repeated here)
  - No sink filling or pit removal is applied to the DEM. Residual SRTM pits
    produce per-scene Max_Depth_m values of 26-37 m, which are artefacts. Use
    Avg_Depth_m; do NOT treat Max_Depth_m as a physical depth.
  - DEM nodata pixels are excluded from the area, the boundary and the volume.
  - Water extent and water depth are NOT independent: the mask determines both
    the area and, through its boundary pixels, the surface elevation from which
    depth is derived. A systematic mask error propagates into the volume twice.
  - There is no independent (bathymetric or GNSS-levelled) validation of the
    depths or volumes produced by this routine.

USAGE
  python scripts/water_volume_from_dem.py \
      --masks   /path/to/water_masks/Chott_Djerid \
      --dem     /path/to/srtm_30m_clipped.tif \
      --out     outputs/Chott_Djerid_water_volumes.csv

  Masks are single-band GeoTIFFs, nonzero = water, one file per date, named so
  that the last 10 characters of the stem are the ISO date, e.g.
  "clipped_masked_-0.1_2019-01-15.tif". Mask and DEM must share the same CRS,
  grid and extent (the masks were reprojected to the DEM grid beforehand).

REQUIREMENTS
  pip install numpy scipy rasterio pandas
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import rasterio
from scipy.ndimage import label, binary_erosion

MIN_PIXELS = 10          # components smaller than this are discarded as noise


def volume_for_one_date(mask_path, dem_data, dem_nodata, pixel_area):
    """Return the per-date summary for a single water mask."""
    with rasterio.open(mask_path) as src:
        mask = src.read(1) > 0

    if mask.shape != dem_data.shape:
        raise ValueError(f"{mask_path.name}: mask and DEM grids differ "
                         f"({mask.shape} vs {dem_data.shape}). Reproject the "
                         f"mask onto the DEM grid first.")

    labeled, n = label(mask)                 # default: 4-connectivity
    total_volume = 0.0
    total_area_m2 = 0.0
    depth_sum = 0.0
    depth_pixels = 0
    max_depth = 0.0
    surface_elevs = []
    bodies = 0

    for obj_id in range(1, n + 1):
        obj_mask = labeled == obj_id
        if np.sum(obj_mask) < MIN_PIXELS:
            continue

        # Boundary for surface elevation estimate
        eroded = binary_erosion(obj_mask)
        boundary = obj_mask ^ eroded

        elev_boundary = dem_data[boundary & (dem_data != dem_nodata)]
        if elev_boundary.size == 0:
            continue

        surface_elev = float(np.mean(elev_boundary))
        surface_elevs.append(surface_elev)

        valid_area = obj_mask & (dem_data != dem_nodata)
        pixel_count = int(np.sum(valid_area))
        surface_area_m2 = pixel_count * pixel_area

        depth = np.where(valid_area, surface_elev - dem_data, 0.0)
        depth = np.where(depth > 0, depth, 0.0)
        volume = float(np.sum(depth * pixel_area))

        total_volume += volume
        total_area_m2 += surface_area_m2
        d = depth[valid_area & (depth > 0)]
        if d.size:
            depth_sum += float(d.sum())
            depth_pixels += int(d.size)
            max_depth = max(max_depth, float(d.max()))
        bodies += 1

    return {
        "Date": mask_path.stem,
        "Water_Volume_m3": round(total_volume, 2),
        "Water_Bodies_Count": bodies,
        "Avg_Depth_m": round(depth_sum / depth_pixels, 3) if depth_pixels else 0.0,
        "Max_Depth_m": round(max_depth, 3),
        "Surface_Elevation_m": round(float(np.mean(surface_elevs)), 2) if surface_elevs else np.nan,
        "Total_Surface_Area_km2": round(total_area_m2 / 1e6, 4),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--masks", required=True, help="directory of water-mask GeoTIFFs")
    ap.add_argument("--dem", required=True, help="DEM GeoTIFF on the same grid as the masks")
    ap.add_argument("--out", required=True, help="output CSV path")
    ap.add_argument("--pattern", default="*.tif", help="mask filename glob (default *.tif)")
    args = ap.parse_args()

    with rasterio.open(args.dem) as src:
        dem_data = src.read(1).astype("float64")
        dem_nodata = src.nodata if src.nodata is not None else -32768
        pixel_area = abs(src.transform.a * src.transform.e)   # m2 per pixel

    masks = sorted(Path(args.masks).glob(args.pattern))
    if not masks:
        raise SystemExit(f"no masks matching {args.pattern} in {args.masks}")

    rows = []
    for i, m in enumerate(masks, 1):
        rows.append(volume_for_one_date(m, dem_data, dem_nodata, pixel_area))
        print(f"[{i}/{len(masks)}] {m.name} -> "
              f"{rows[-1]['Water_Volume_m3']:.0f} m3, "
              f"{rows[-1]['Total_Surface_Area_km2']:.2f} km2")

    out = Path(args.out); out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
