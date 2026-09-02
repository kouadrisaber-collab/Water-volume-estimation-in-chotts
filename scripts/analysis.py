# -*- coding: utf-8 -*-
"""Round-3 revision analyses: trend tests, climate-water correlations, Li sensitivity."""
# --- repository paths (resolve relative to this file, not the shell's cwd) ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTDIR = ROOT / "outputs"
OUTDIR.mkdir(exist_ok=True)
# ---------------------------------------------------------------------------

import glob, os, re, json
import numpy as np, pandas as pd
from scipy import stats
import pymannkendall as mk

OUT = {}
FILES = {
 'Djerid':  DATA / "Chott_Djerid_water_volumes.csv",
 'Gharsa':  DATA / "Chott_Gharsa_water_volumes.csv",
 'Melghigh':DATA / "Chott_Melghigh_water_volumes.csv",
 'Merouane':DATA / "Chott_Merouan_water_volumes.csv"}

# ---- 0. thresholds actually used per sensor -------------------------------
thr = {}
for c,f in FILES.items():
    d = pd.read_csv(f)
    thr[c] = sorted(set(re.findall(r'(-?\d+\.\d+)', ' '.join(d['Date'].astype(str).str.rsplit('_',n=1).str[0]))))
OUT['thresholds_used'] = thr

# ---- 1. load & monthly-aggregate water series -----------------------------
recs = []
for c,f in FILES.items():
    d = pd.read_csv(f)
    d['date'] = pd.to_datetime(d['Date'].str[-10:], errors='coerce')
    d = d.dropna(subset=['date'])
    d['chott'] = c
    recs.append(d[['chott','date','Water_Volume_m3','Avg_Depth_m','Max_Depth_m',
                   'Surface_Elevation_m','Total_Surface_Area_km2','Water_Bodies_Count']])
w = pd.concat(recs)
w['year'] = w.date.dt.year; w['month'] = w.date.dt.month
wm = (w.groupby(['chott','year','month'])
        .agg(area=('Total_Surface_Area_km2','mean'), vol=('Water_Volume_m3','mean'),
             depth=('Avg_Depth_m','mean'), nscenes=('date','size')).reset_index())
OUT['n_scenes'] = {c:int((w.chott==c).sum()) for c in FILES}
OUT['scene_span'] = {c: [str(w[w.chott==c].date.min().date()), str(w[w.chott==c].date.max().date())] for c in FILES}

# ---- 2. climate merge -----------------------------------------------------
cl = pd.read_csv(DATA / "chotts_climate_enriched.csv")
cl['chott'] = cl['chott'].str.capitalize().replace({'Merouan':'Merouane'})
cl = cl[['chott','year','month','rain_mm','temp_C','pet_mm','rh_pct','lst_C']]
m = wm.merge(cl, on=['chott','year','month'], how='inner').sort_values(['chott','year','month'])

# ---- 3. correlations with lagged rainfall ---------------------------------
corr_rows = []
for c, g in m.groupby('chott'):
    g = g.reset_index(drop=True)
    for lag in range(0, 7):
        # cumulative antecedent rainfall over `lag` preceding months (lag=0 -> same month)
        rain_cum = g['rain_mm'].rolling(lag+1, min_periods=lag+1).sum()
        for var in ['area','vol']:
            ok = rain_cum.notna() & g[var].notna()
            if ok.sum() < 12: continue
            r,pr = stats.pearsonr(rain_cum[ok], g[var][ok])
            rho,ps = stats.spearmanr(rain_cum[ok], g[var][ok])
            corr_rows.append(dict(chott=c, var=var, predictor='rain_cum', lag_months=lag, n=int(ok.sum()),
                                  pearson_r=round(r,3), p_pearson=round(pr,4),
                                  spearman_rho=round(rho,3), p_spearman=round(ps,4)))
    for var in ['area','vol']:
        for pred in ['pet_mm','temp_C','rh_pct']:
            ok = g[pred].notna() & g[var].notna()
            r,pr = stats.pearsonr(g[pred][ok], g[var][ok]); rho,ps = stats.spearmanr(g[pred][ok], g[var][ok])
            corr_rows.append(dict(chott=c, var=var, predictor=pred, lag_months=0, n=int(ok.sum()),
                                  pearson_r=round(r,3), p_pearson=round(pr,4),
                                  spearman_rho=round(rho,3), p_spearman=round(ps,4)))
corr = pd.DataFrame(corr_rows)
corr.to_csv(OUTDIR / "correlations.csv", index=False)

# ---- 4. trend tests -------------------------------------------------------
trend_rows = []
def add_trend(chott, series_name, x, kind='mk'):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if len(x) < 4: return
    res = mk.seasonal_test(x, period=12) if kind=='seasonal' else mk.original_test(x)
    trend_rows.append(dict(chott=chott, series=series_name, test=('Seasonal MK' if kind=='seasonal' else 'MK'),
                           n=len(x), trend=res.trend, tau=round(res.Tau,3), p=round(res.p,4),
                           sens_slope=round(res.slope,4), significant_05=bool(res.p<0.05)))
for c, g in m.groupby('chott'):
    g = g.sort_values(['year','month'])
    add_trend(c, 'monthly water area (km2)', g['area'], 'seasonal')
    add_trend(c, 'monthly water volume (m3)', g['vol'], 'seasonal')
    for mo,lab in [(2,'wet-season (Feb) water area (km2)'),(6,'dry-season (Jun) water area (km2)')]:
        s = g[g.month==mo].sort_values('year')['area']
        add_trend(c, lab, s)
    ann = g.groupby('year')['area'].mean()
    add_trend(c, 'annual mean water area (km2)', ann)
    annv = g.groupby('year')['vol'].mean()
    add_trend(c, 'annual mean water volume (m3)', annv)
trend = pd.DataFrame(trend_rows)
trend.to_csv(OUTDIR / "trends.csv", index=False)

# ---- 5. lithium concentration sensitivity --------------------------------
# mean annual volume per chott per season, from the manuscript's Feb/Jun analysis
seas = {}
for c, g in m.groupby('chott'):
    seas[c] = dict(wet_mean_vol=float(g[g.month==2]['vol'].mean()),
                   dry_mean_vol=float(g[g.month==6]['vol'].mean()),
                   all_mean_vol=float(g['vol'].mean()))
OUT['seasonal_mean_volumes_m3'] = seas
base = {'Djerid':44.1,'Gharsa':44.1,'Melghigh':66.0,'Merouane':66.0}
mults = [0.5,0.75,1.0,1.5,2.0]
sens_rows=[]
for c in FILES:
    for mu in mults:
        C = base[c]*mu   # mg/L
        # tonnes = V(m3) * 1000 L/m3 * C(mg/L) * 1e-9 t/mg
        for lab,V in [('wet(Feb)',seas[c]['wet_mean_vol']),('dry(Jun)',seas[c]['dry_mean_vol'])]:
            sens_rows.append(dict(chott=c, season=lab, C_mgL=round(C,1), multiplier=mu,
                                  mean_volume_m3=round(V), Li_tonnes=round(V*1000*C*1e-9)))
sens = pd.DataFrame(sens_rows); sens.to_csv(OUTDIR / "li_sensitivity.csv", index=False)

# ---- 6. boundary-pixel / depth descriptive stats --------------------------
bp = w.groupby('chott').agg(mean_avg_depth=('Avg_Depth_m','mean'), sd_avg_depth=('Avg_Depth_m','std'),
        mean_max_depth=('Max_Depth_m','mean'), mean_surf_elev=('Surface_Elevation_m','mean'),
        sd_surf_elev=('Surface_Elevation_m','std'), mean_bodies=('Water_Bodies_Count','mean')).round(3)
bp.to_csv(OUTDIR / "depth_stats.csv")
OUT['depth_stats'] = bp.reset_index().to_dict('records')

json.dump(OUT, open(OUTDIR / "summary.json",'w'), indent=1, default=str)
print(json.dumps(OUT, indent=1, default=str)[:2500])
