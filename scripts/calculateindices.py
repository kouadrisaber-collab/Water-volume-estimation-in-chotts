# -*- coding: utf-8 -*-
"""CalculateIndices.ipynb"""

import ee
import pandas as pd
import os

# ============================
# Authenticate Earth Engine
# ============================
ee.Authenticate()   # Follow link, copy-paste token
ee.Initialize(project='ee-kouadrisaber8')

from google.colab import drive

# Classification threshold applied to MNDSI = (SWIR1 - SWIR2)/(SWIR1 + SWIR2).
MNDSI_THRESHOLD = 0.12

drive.mount('/content/drive', force_remount=True)

# ---------------------------
# Regions of Interest (Chotts)
# ---------------------------
chotts = {
    'Chott_Melghigh': ee.Geometry.Rectangle([5.825829, 34.089807, 6.996115, 34.913205]),
    'Chott_Merouane': ee.Geometry.Rectangle([5.927028, 33.853596, 6.500225, 34.212210]),
    'Chott_Djerid': ee.Geometry.Rectangle([7.724641, 33.290255, 9.464564, 34.065671]),
    'Chott_Gharsa': ee.Geometry.Rectangle([7.360559, 33.915131, 8.235854, 34.258467])
}

# ============================
# 1. Landsat 8 Processing
# ============================
print("Starting Landsat 8 Processing...")
l8_start_date = '2019-01-01'
l8_end_date   = '2024-12-31'
l8_images_per_month = 1
l8_min_coverage = 0.95

def mask_clouds_landsat(image):
    qa = image.select('QA_PIXEL')
    cloud_free = qa.bitwiseAnd(1 << 5).eq(0)
    shadow_free = qa.bitwiseAnd(1 << 3).eq(0)
    return image.updateMask(cloud_free.And(shadow_free))

def add_indices_landsat(image):
    swir1 = image.select('SR_B6')
    swir2 = image.select('SR_B7')
    # Calculate MNDSI
    mndsi = swir1.subtract(swir2).divide(swir1.add(swir2)).rename('MNDSI')
    return image.addBands([mndsi])

def add_coverage_property_landsat(region_geom):
    def wrap(img):
        footprint = ee.Image(img).geometry()
        intersection = footprint.intersection(region_geom, ee.ErrorMargin(1))
        cov = intersection.area().divide(region_geom.area())
        return img.set('coverage', cov)
    return wrap

l8_summary = []

for chott_name, region_geom in chotts.items():
    print(f"\n--- Processing Landsat 8 for {chott_name} ---")
    
    collection = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
        .filterDate(l8_start_date, l8_end_date)
        .filterBounds(region_geom)
        .map(add_coverage_property_landsat(region_geom))
        .filter(ee.Filter.gte('coverage', l8_min_coverage))
        .map(mask_clouds_landsat)
        .map(add_indices_landsat)
        .sort('system:time_start')
    )

    col_size = collection.size().getInfo()
    print(f"Total images after filtering for {chott_name}: {col_size}")

    if col_size > 0:
        images_list = collection.toList(col_size)
        selected_images = []

        for i in range(col_size):
            img = ee.Image(images_list.get(i))
            date_str = ee.Date(img.get('system:time_start')).format('YYYY-MM-dd').getInfo()
            year = ee.Date(img.get('system:time_start')).get('year').getInfo()
            month = ee.Date(img.get('system:time_start')).get('month').getInfo()
            selected_images.append((year, month, date_str, img))

        df_imgs = pd.DataFrame(selected_images, columns=['year','month','date','image'])
        df_imgs['ym'] = df_imgs['year'].astype(str) + "-" + df_imgs['month'].astype(str).str.zfill(2)

        selected_rows = []
        for ym, group in df_imgs.groupby('ym'):
            pick = group.sort_values('date').head(l8_images_per_month)
            for _, row in pick.iterrows():
                selected_rows.append((row['date'], row['image']))

        print(f"Selected {len(selected_rows)} images.")

        for date_str, img in selected_rows:
            image = ee.Image(img).clip(region_geom)
            
            mask = image.select('MNDSI').gt(MNDSI_THRESHOLD)
            pixel_area = ee.Image.pixelArea().updateMask(mask)
            stats = pixel_area.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=region_geom,
                scale=30,
                maxPixels=1e13,
                bestEffort=True
            )
            try:
                area_m2 = stats.getInfo().get('area', 0)
                if area_m2 is None:
                    area_m2 = 0
            except:
                area_m2 = 0
                
            area_km2 = round(area_m2 / 1e6, 4)
            l8_summary.append({
                'Date': date_str,
                'Chott': chott_name,
                'MNDSI_km2': area_km2
            })
            print(f"{date_str} -> MNDSI_km2: {area_km2}")
    else:
        print(f"⚠️ No images passed the filters for {chott_name}.")

# Save Landsat 8 results as a pivot table
if l8_summary:
    df_l8 = pd.DataFrame(l8_summary)
    # Pivot so each Chott is a column and Dates are rows
    df_l8_pivot = df_l8.pivot_table(index='Date', columns='Chott', values='MNDSI_km2').reset_index()
    l8_output_path = '/content/drive/MyDrive/AllChotts_Landsat8_MNDSI.xlsx'
    df_l8_pivot.to_excel(l8_output_path, index=False)
    print(f"\n✅ Landsat 8 Results saved to {l8_output_path}")


# ============================
# 2. Sentinel-2 Processing
# ============================
print("\nStarting Sentinel-2 Processing...")
s2_start_date = '2020-01-01'
s2_end_date   = '2020-06-30'
s2_images_per_month = 2
s2_min_coverage = 0.5

def mask_clouds_s2(image):
    qa = image.select('QA60')
    cloud_bit_mask  = int('1000000000000', 2)
    cirrus_bit_mask = int('10000000000', 2)
    mask = qa.bitwiseAnd(cloud_bit_mask).eq(0).And(
           qa.bitwiseAnd(cirrus_bit_mask).eq(0))
    return image.updateMask(mask).divide(10000)

def add_indices_s2(image):
    swir1 = image.select('B11')
    swir2 = image.select('B12')
    mndsi = swir1.subtract(swir2).divide(swir1.add(swir2)).rename('MNDSI')
    return image.addBands([mndsi])

def add_coverage_property_s2(region_geom):
    def wrap(img):
        footprint = ee.Image(img).geometry()
        intersection = footprint.intersection(region_geom, ee.ErrorMargin(1))
        cov = intersection.area().divide(region_geom.area())
        return img.set('coverage', cov)
    return wrap

s2_summary = []

for chott_name, region_geom in chotts.items():
    print(f"\n--- Processing Sentinel-2 for {chott_name} ---")
    
    collection = (
        ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
        .filterDate(s2_start_date, s2_end_date)
        .filterBounds(region_geom)
        .map(add_coverage_property_s2(region_geom))
        .filter(ee.Filter.gte('coverage', s2_min_coverage))
        .map(mask_clouds_s2)
        .map(add_indices_s2)
        .sort('system:time_start')
    )

    col_size = collection.size().getInfo()
    print(f"Total images after filtering for {chott_name}: {col_size}")

    if col_size > 0:
        images_list = collection.toList(col_size)
        selected_images = []

        for i in range(col_size):
            img_obj = images_list.get(i)
            if img_obj is None:
                continue
            img = ee.Image(img_obj)
            date_val = img.get('system:time_start')
            if date_val is None:
                continue
            date_str = ee.Date(date_val).format('YYYY-MM-dd').getInfo()
            year = ee.Date(date_val).get('year').getInfo()
            month = ee.Date(date_val).get('month').getInfo()
            selected_images.append((year, month, date_str, img))

        df_imgs = pd.DataFrame(selected_images, columns=['year','month','date','image'])
        df_imgs['ym'] = df_imgs['year'].astype(str) + "-" + df_imgs['month'].astype(str).str.zfill(2)

        selected_rows = []
        for ym, group in df_imgs.groupby('ym'):
            pick = group.sort_values('date').head(s2_images_per_month)
            for _, row in pick.iterrows():
                selected_rows.append((row['date'], row['image']))

        print(f"Selected {len(selected_rows)} images.")

        for date_str, img in selected_rows:
            image = ee.Image(img).clip(region_geom)
            mask = image.select('MNDSI').gt(MNDSI_THRESHOLD)
            pixel_area = ee.Image.pixelArea().updateMask(mask)
            stats = pixel_area.reduceRegion(
                reducer=ee.Reducer.sum(),
                geometry=region_geom,
                scale=20,
                maxPixels=1e13,
                bestEffort=True
            )
            try:
                area_m2 = stats.getInfo().get('area', 0)
                if area_m2 is None:
                    area_m2 = 0
            except:
                area_m2 = 0
                
            area_km2 = round(area_m2 / 1e6, 4)
            s2_summary.append({
                'Date': date_str,
                'Chott': chott_name,
                'MNDSI_km2': area_km2
            })
            print(f"{date_str} -> MNDSI_km2: {area_km2}")
    else:
        print(f"⚠️ No images passed the filters for {chott_name}.")

# Save Sentinel-2 results as a pivot table
if s2_summary:
    df_s2 = pd.DataFrame(s2_summary)
    df_s2_pivot = df_s2.pivot_table(index='Date', columns='Chott', values='MNDSI_km2').reset_index()
    s2_output_path = '/content/drive/MyDrive/AllChotts_Sentinel2_MNDSI.xlsx'
    df_s2_pivot.to_excel(s2_output_path, index=False)
    print(f"\n✅ Sentinel-2 Results saved to {s2_output_path}")
