
# --- repository paths (resolve relative to this file, not the shell's cwd) ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTDIR = ROOT / "outputs"
OUTDIR.mkdir(exist_ok=True)
# ---------------------------------------------------------------------------

import pandas as pd, numpy as np, openpyxl
from scipy import stats
import pymannkendall as mk

d = pd.read_excel(DATA / "ChottDjerid_Sentinel2_MNDSI.xlsx")
d = d.rename(columns={'MNDSI_km2':'sas'})[['Date','sas']]; d['chott']='Djerid'
o = pd.read_excel(DATA / "Marouane, melghigh, gharssa Chotts_Landsat8_MNDSI.xlsx")
o = o.melt(id_vars='Date', var_name='chott', value_name='sas').dropna()
o['chott'] = o.chott.str.replace('Chott_','')
sas = pd.concat([d,o]); sas['date']=pd.to_datetime(sas.Date)
sas['year']=sas.date.dt.year; sas['month']=sas.date.dt.month
sm = sas.groupby(['chott','year','month'])['sas'].mean().reset_index()

FILES={'Djerid':'Chott_Djerid','Gharsa':'Chott_Gharsa','Melghigh':'Chott_Melghigh','Merouane':'Chott_Merouan'}
recs=[]
for c,f in FILES.items():
    x=pd.read_csv(DATA / f"{f}_water_volumes.csv")
    x['date']=pd.to_datetime(x['Date'].str[-10:],errors='coerce'); x=x.dropna(subset=['date'])
    x['chott']=c; x['year']=x.date.dt.year; x['month']=x.date.dt.month
    recs.append(x)
w=pd.concat(recs).groupby(['chott','year','month']).agg(area=('Total_Surface_Area_km2','mean'),vol=('Water_Volume_m3','mean')).reset_index()
j=w.merge(sm,on=['chott','year','month'])

rows=[]
for c,g in j.groupby('chott'):
    g=g.sort_values(['year','month'])
    r,pr=stats.pearsonr(g.area,g.sas); rho,ps=stats.spearmanr(g.area,g.sas)
    t=mk.seasonal_test(g.sas.values,period=12); ta=mk.original_test(g.groupby('year')['sas'].mean().values)
    rows.append(dict(chott=c,n=len(g),water_vs_SAS_pearson=round(r,3),p_pearson=round(pr,4),
        water_vs_SAS_spearman=round(rho,3),p_spearman=round(ps,4),
        SAS_seasonalMK_trend=t.trend,SAS_tau=round(t.Tau,3),SAS_p=round(t.p,4),SAS_sens_slope=round(t.slope,3),
        SAS_annualMK_trend=ta.trend,SAS_annual_p=round(ta.p,4),
        SAS_mean_km2=round(g.sas.mean(),1)))
out=pd.DataFrame(rows); out.to_csv(OUTDIR / "salinity_trends_and_coupling.csv",index=False)
print(out.T.to_string())
