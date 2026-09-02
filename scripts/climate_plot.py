# =====================================================================
#  Build the supplementary climate figure (4-panel climograph + water).
#  Input : chotts_climate_2019_2024.csv  (from climate_extract_gee.py)
#  Output: Fig_S1_climate.png / .pdf
#  Requires: pandas, matplotlib, numpy
# =====================================================================
# --- repository paths (resolve relative to this file, not the shell's cwd) ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTDIR = ROOT / "outputs"
OUTDIR.mkdir(exist_ok=True)
# ---------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

CSV = DATA / "chotts_climate_2019_2024.csv"

# ---------------------------------------------------------------------
# Water-surface area (km2) from Table 3 of the manuscript.
# Plotted at Feb (wet, month=2) and June (dry, month=6) of each year.
# ---------------------------------------------------------------------
water = {
 'Djerid':   {'wet': {2019:735.33,2020:3451.53,2021:705.02,2022:1098.67,2023:309.95,2024:302.72},
              'dry': {2019:1501.9,2020:345.74,2021:2913.69,2022:296.53,2023:450.72,2024:149.74}},
 'Gharsa':   {'wet': {2019:188.73,2020:674.63,2021:117.29,2022:265.34,2023:24.03,2024:41.82},
              'dry': {2019:807.74,2020:177.10,2021:138.99,2022:46.84,2023:146.70,2024:114.85}},
 'Melghigh': {'wet': {2019:540.15,2020:502.39,2021:491.02,2022:143.45,2023:63.05,2024:60.88},
              'dry': {2019:392.27,2020:378.79,2021:1005.45,2022:64.35,2023:264.25,2024:286.20}},
 'Merouane': {'wet': {2019:353.03,2020:421.30,2021:393.28,2022:214.92,2023:230.35,2024:207.38},
              'dry': {2019:474.25,2020:449.12,2021:474.33,2022:288.62,2023:276.69,2024:259.10}},
}
ORDER = ['Djerid', 'Gharsa', 'Melghigh', 'Merouane']

# ---------------------------------------------------------------------
# Load + unit conversions
# ---------------------------------------------------------------------
df = pd.read_csv(CSV)
df['chott'] = df['chott'].str.strip().str.capitalize()   # normalise e.g. 'melghigh' -> 'Melghigh'
df['t'] = pd.to_datetime(df['date'] + '-01')
df['temp_C'] = df['t2m_K'] - 273.15
df['pet_mm'] = df['pet_raw'] * 0.1
df['aet_mm'] = df['aet_raw'] * 0.1
df['lst_C']  = df['lst_raw'] * 0.02 - 273.15
# Relative humidity (%) from air temperature and dewpoint via the Magnus formula.
# RH = 100 * e_s(Td) / e_s(T), with e_s(T) = exp(17.625*T / (243.04+T)), T in degC.
if 'd2m_K' in df.columns:
    Tc  = df['t2m_K'] - 273.15
    Tdc = df['d2m_K'] - 273.15
    df['rh_pct'] = (100.0 * np.exp(17.625 * Tdc / (243.04 + Tdc))
                          / np.exp(17.625 * Tc  / (243.04 + Tc))).clip(0, 100)
df = df.sort_values(['chott', 't'])
# Save an enriched table (with converted units + RH) alongside the figure.
df.to_csv(OUTDIR / "chotts_climate_enriched.csv", index=False)

def offset_axis(ax, offset):
    """Create a twin y-axis with its spine pushed out to the right."""
    a = ax.twinx()
    a.spines['right'].set_position(('outward', offset))
    return a

fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True)
fig.subplots_adjust(left=0.11, right=0.88, hspace=0.22, wspace=0.95)

def left_offset_axis(ax, offset):
    """Create a twin y-axis pushed out to the LEFT of the plot."""
    a = ax.twinx()
    a.spines['left'].set_position(('outward', offset))
    a.yaxis.set_label_position('left')
    a.yaxis.set_ticks_position('left')
    a.spines['right'].set_visible(False)
    return a

for ax1, name in zip(axes.flat, ORDER):
    d = df[df['chott'] == name]

    # --- LEFT (inner) axis : rainfall bars on their own scale (mm) ---
    ax1.bar(d['t'], d['rain_mm'], width=22, color='#4a90d9', alpha=0.9)
    ax1.set_ylabel('Rainfall (mm)', color='#2b6cb0')
    ax1.tick_params(axis='y', labelcolor='#2b6cb0')
    ax1.set_title('Chott ' + name, fontweight='bold')
    ax1.margins(x=0.01)

    # --- LEFT (outer) axis : PET line on its own scale (mm) ---
    axp = left_offset_axis(ax1, 46)
    axp.plot(d['t'], d['pet_mm'], color='#d1495b', lw=1.4, ls='--')
    axp.set_ylabel('PET (mm)', color='#d1495b')
    axp.tick_params(axis='y', labelcolor='#d1495b')

    # --- RIGHT (inner) axis : air temperature + LST (deg C) ---
    ax2 = ax1.twinx()
    ax2.plot(d['t'], d['temp_C'], color='#e08e0b', lw=1.6, label='Air temp (ERA5-Land)')
    ax2.plot(d['t'], d['lst_C'], color='#7a4019', lw=1.0, ls=':', label='LST day (MODIS)')
    ax2.set_ylabel('Temperature (°C)', color='#e08e0b')
    ax2.tick_params(axis='y', labelcolor='#e08e0b')

    # --- FAR-RIGHT offset axis : water surface area (km2) ---
    ax3 = offset_axis(ax1, 48)
    wy = sorted(water[name]['wet'])
    ax3.scatter([pd.Timestamp(y, 2, 1) for y in wy],
                [water[name]['wet'][y] for y in wy],
                marker='o', s=42, facecolor='#2e8b57', edgecolor='k',
                zorder=5, label='Water area, wet (Feb)')
    dy = sorted(water[name]['dry'])
    ax3.scatter([pd.Timestamp(y, 6, 1) for y in dy],
                [water[name]['dry'][y] for y in dy],
                marker='s', s=42, facecolor='#7b3fa0', edgecolor='k',
                zorder=5, label='Water area, dry (Jun)')
    ax3.set_ylabel('Water area (km²)', color='#2e8b57')
    ax3.tick_params(axis='y', labelcolor='#2e8b57')

    ax1.xaxis.set_major_locator(mdates.YearLocator())
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

# Single shared legend built from fixed proxy handles (clean and order-stable)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_items = [
    Patch(facecolor='#4a90d9', alpha=0.85, label='Rainfall (CHIRPS)'),
    Line2D([0], [0], color='#d1495b', lw=1.4, ls='--', label='PET (TerraClimate)'),
    Line2D([0], [0], color='#e08e0b', lw=1.6, label='Air temp (ERA5-Land)'),
    Line2D([0], [0], color='#7a4019', lw=1.0, ls=':', label='LST day (MODIS)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#2e8b57',
           markeredgecolor='k', markersize=8, label='Water area, wet (Feb)'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor='#7b3fa0',
           markeredgecolor='k', markersize=8, label='Water area, dry (Jun)'),
]
fig.legend(handles=legend_items, loc='lower center', ncol=6, frameon=False,
           bbox_to_anchor=(0.5, -0.02))

fig.suptitle('Supplementary climate context (2019–2024): rainfall, evaporative demand, '
             'temperature and water-surface area', fontweight='bold', y=0.98)
fig.savefig(OUTDIR / "Fig09_climate_context.png", dpi=300, bbox_inches='tight')
fig.savefig(OUTDIR / "Fig09_climate_context.pdf", bbox_inches='tight')
print('Saved Fig09_climate_context.png/.pdf to', OUTDIR)
