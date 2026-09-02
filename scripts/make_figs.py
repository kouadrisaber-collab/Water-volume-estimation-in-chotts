
# --- repository paths (resolve relative to this file, not the shell's cwd) ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTDIR = ROOT / "outputs"
OUTDIR.mkdir(exist_ok=True)
# ---------------------------------------------------------------------------

import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9})
CH=['Djerid','Gharsa','Melghigh','Merouane']; COL=dict(zip(CH,['#1f77b4','#d62728','#2ca02c','#9467bd']))

c=pd.read_csv(OUTDIR / "correlations.csv")
fig,axes=plt.subplots(1,4,figsize=(11,3.0),sharey=True)
for ax,ch in zip(axes,CH):
    g=c[(c.chott==ch)&(c['var']=='area')&(c.predictor=='rain_cum')].sort_values('lag_months')
    ax.axhline(0,color='0.6',lw=.8)
    ax.plot(g.lag_months,g.pearson_r,'-o',color=COL[ch],ms=5,label='Pearson r')
    ax.plot(g.lag_months,g.spearman_rho,'--s',color=COL[ch],ms=4,mfc='white',label='Spearman ρ')
    for _,r in g.iterrows():
        if r.p_pearson<0.05:
            ax.annotate('*',(r.lag_months,r.pearson_r),textcoords='offset points',xytext=(0,7),
                        ha='center',fontsize=12,color=COL[ch])
    ax.set_title(f'Chott {ch}',fontsize=10); ax.set_xlabel('Antecedent rainfall window (months)')
    ax.set_ylim(-0.25,0.75); ax.grid(alpha=.25,lw=.5)
axes[0].set_ylabel('Correlation with monthly\nwater area')
axes[0].legend(fontsize=7.5,frameon=False,loc='upper left')
fig.suptitle('Correlation between cumulative antecedent rainfall (CHIRPS) and monthly water area, 2019–2024  (* p < 0.05)',
             fontsize=9.5,y=1.02)
fig.tight_layout(); fig.savefig(OUTDIR / "Fig07_rainfall_lag_correlation.png",dpi=400,bbox_inches='tight')

s=pd.read_csv(OUTDIR / "li_sensitivity.csv")
fig,axes=plt.subplots(1,4,figsize=(11,3.0),sharey=False)
for ax,ch in zip(axes,CH):
    g=s[s.chott==ch]
    for seas,mk_,ls in [('wet(Feb)','o','-'),('dry(Jun)','s','--')]:
        gg=g[g.season==seas].sort_values('C_mgL')
        ax.plot(gg.C_mgL,gg.Li_tonnes/1000,ls,marker=mk_,color=COL[ch],ms=5,
                mfc='white' if seas=='dry(Jun)' else COL[ch],label='Wet (Feb)' if seas=='wet(Feb)' else 'Dry (Jun)')
    base=44.1 if ch in('Djerid','Gharsa') else 66.0
    ax.axvline(base,color='0.4',ls=':',lw=1)
    ax.annotate(f'literature\n{base} mg L$^{{-1}}$',(base,ax.get_ylim()[1]*0.92),fontsize=7,ha='center',color='0.3')
    ax.set_title(f'Chott {ch}',fontsize=10); ax.set_xlabel('Assumed Li concentration (mg L$^{-1}$)')
    ax.grid(alpha=.25,lw=.5)
axes[0].set_ylabel('Lithium potential (10³ t)'); axes[0].legend(fontsize=7.5,frameon=False)
fig.suptitle('Sensitivity of estimated lithium potential to the assumed brine lithium concentration (0.5×–2.0× the literature value)',
             fontsize=9.5,y=1.02)
fig.tight_layout(); fig.savefig(OUTDIR / "Fig08_lithium_concentration_sensitivity.png",dpi=400,bbox_inches='tight')
print('figures written')
