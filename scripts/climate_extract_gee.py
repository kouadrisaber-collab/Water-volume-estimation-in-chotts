# =====================================================================
#  Supplementary climate data extraction for the Algero-Tunisian chotts
#  Run this in your existing Google Colab / GEE Python environment.
#  Output: chotts_climate_2019_2024.csv  (one row per chott per month)
# =====================================================================
import ee
ee.Authenticate()                          # opens the Google sign-in flow (run once per session)
ee.Initialize(project='your-project-id')   # <-- REPLACE with your GEE Cloud project ID
#   The project ID is an identifier (see https://code.earthengine.google.com, project
#   dropdown), not a secret. If ee.Authenticate() says credentials already exist, it
#   simply skips the prompt.

# ---------------------------------------------------------------------
# 1. STUDY PERIOD : 72 months, Jan 2019 -> Dec 2024
# ---------------------------------------------------------------------
START  = ee.Date('2019-01-01')
MONTHS = ee.List.sequence(0, 71)

# ---------------------------------------------------------------------
# 2. CHOTT GEOMETRIES  ***  REPLACE THESE WITH YOUR OWN POLYGONS  ***
#    Use the SAME geometries you already use for the imagery analysis,
#    e.g.  'Djerid': ee.FeatureCollection('users/you/chott_djerid').geometry()
#    The rectangles below are only rough placeholders.
# ---------------------------------------------------------------------
chotts = {
    'melghigh': ee.Geometry.Rectangle([5.825829, 34.089807, 6.996115, 34.913205]),
    'Merouane': ee.Geometry.Rectangle([5.927028, 33.853596, 6.500225, 34.212210]),  
    'Djerid': ee.Geometry.Rectangle([7.724641, 33.290255, 9.464564, 34.065671]),  
    'Gharsa': ee.Geometry.Rectangle([7.360559, 33.915131, 8.235854, 34.258467]),  
}

# ---------------------------------------------------------------------
# 3. SOURCE COLLECTIONS  (native scales kept to limit computation)
# ---------------------------------------------------------------------
chirps = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY').select('precipitation')      # mm/day
# ERA5-Land carries no relative-humidity band; we export temperature AND dewpoint
# and derive relative humidity (Magnus formula) in the plotting script.
era5   = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') \
            .select(['temperature_2m', 'dewpoint_temperature_2m'])               # K
tclim  = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').select(['pet', 'aet'])   # *0.1 = mm

# MODIS LST is an 8-day product that is ALREADY cloud-screened (clear-sky only),
# so no image cloud mask is required here. The only optional quality step is to
# drop low-quality retrievals using the QC_Day band (mandatory QA flag 0 = good,
# 1 = other quality). Set USE_LST_QC = False to keep every retrieved pixel.
USE_LST_QC = True
def _mask_lst(img):
    qc = img.select('QC_Day')
    good = qc.bitwiseAnd(3).lte(1)                  # mandatory QA in {0, 1}
    return img.select('LST_Day_1km').updateMask(good)
modis = ee.ImageCollection('MODIS/061/MOD11A2')
modis = modis.map(_mask_lst) if USE_LST_QC else modis.select('LST_Day_1km')      # *0.02 = K

# Raw values are exported; unit conversions are done in the plotting script
# (this avoids null-arithmetic errors when a dataset is missing for a month).
def monthly_for_chott(name, geom):
    def per_month(m):
        d0 = START.advance(m, 'month')
        d1 = d0.advance(1, 'month')

        rain = chirps.filterDate(d0, d1).sum() \
            .reduceRegion(ee.Reducer.mean(), geom, 5566, maxPixels=1e13).get('precipitation')
        e5 = era5.filterDate(d0, d1).mean() \
            .reduceRegion(ee.Reducer.mean(), geom, 11132, maxPixels=1e13)
        t2m = e5.get('temperature_2m')
        d2m = e5.get('dewpoint_temperature_2m')
        tc = tclim.filterDate(d0, d1).mean()
        pet = tc.select('pet').reduceRegion(ee.Reducer.mean(), geom, 4638, maxPixels=1e13).get('pet')
        aet = tc.select('aet').reduceRegion(ee.Reducer.mean(), geom, 4638, maxPixels=1e13).get('aet')
        lst = modis.filterDate(d0, d1).mean() \
            .reduceRegion(ee.Reducer.mean(), geom, 1000, maxPixels=1e13).get('LST_Day_1km')

        return ee.Feature(None, {
            'chott':   name,
            'date':    d0.format('YYYY-MM'),
            'year':    d0.get('year'),
            'month':   d0.get('month'),
            'rain_mm': rain,    # monthly total, mm
            't2m_K':   t2m,     # mean air temperature, Kelvin
            'd2m_K':   d2m,     # mean dewpoint temperature, Kelvin (-> RH)
            'pet_raw': pet,     # TerraClimate PET, *0.1 -> mm
            'aet_raw': aet,     # TerraClimate actual ET, *0.1 -> mm
            'lst_raw': lst,     # MODIS day LST, *0.02 -> K
        })
    return ee.FeatureCollection(MONTHS.map(per_month))

# ---------------------------------------------------------------------
# 4. BUILD THE FULL TABLE  (4 chotts x 72 months = 288 rows)
# ---------------------------------------------------------------------
out = ee.FeatureCollection([monthly_for_chott(n, g) for n, g in chotts.items()]).flatten()

COLS = ['chott', 'date', 'year', 'month',
        'rain_mm', 't2m_K', 'd2m_K', 'pet_raw', 'aet_raw', 'lst_raw']

# ---- OPTION A (recommended): asynchronous export to Google Drive ----
task = ee.batch.Export.table.toDrive(
    collection=out,
    description='chotts_climate_2019_2024',
    fileNamePrefix='chotts_climate_2019_2024',
    fileFormat='CSV',
    selectors=COLS,
)
task.start()
print('Export task started. Check the GEE Tasks tab; the CSV will appear in your Drive.')
print('Then download it next to climate_plot.py and run that script.')

# ---- OPTION B (inline, if the export feels slow): direct to pandas ----
# import pandas as pd
# feats = out.getInfo()['features']
# df = pd.DataFrame([f['properties'] for f in feats])
# df.to_csv('chotts_climate_2019_2024.csv', index=False)
# print(df.head())
