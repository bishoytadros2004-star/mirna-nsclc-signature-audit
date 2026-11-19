# Step 8: Figure 1 (published sets vs random sets) from the saved results and nulls.
import os, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES=os.environ.get("RESULTS_DIR","../results")
d=pd.read_csv(os.path.join(RES,"audit_results_by_signature.csv"))
e=pd.read_csv(os.path.join(RES,"external_validation_GSE16025_LUSC.csv")); e=e[e.k>0]
order=["Li2017","Gallach2017","MacDonagh2021","Yu2008","Raponi2009","Hu2010","Guo2016","Zhang2019","Liu2017","Zhou2020","Xue2020","Yang2019","Effusion2012"]
panels=[("TCGA-LUAD",d[d.cohort=="LUAD"],"LUAD"),("TCGA-LUSC",d[d.cohort=="LUSC"],"LUSC"),("GSE16025 LUSC (n=61)",e,None)]
fig,axs=plt.subplots(1,3,figsize=(13,5.2),sharey=True)
for ax,(title,df,coh) in zip(axs,panels):
    df=df.set_index("set").reindex([o for o in order if o in set(df.set)]); y=np.arange(len(df))[::-1]
    for yi,(nm,r) in zip(y,df.iterrows()):
        if coh:
            nulls=json.load(open(f"result_{coh}.json"))["null"][str(int(r.k))]["only"]
            lo,hi=np.nanpercentile(nulls,[2.5,97.5]); ax.plot([lo,hi],[yi,yi],color="0.8",lw=7,solid_capstyle="butt")
        ax.plot(r.null_med,yi,"|",color="0.4",ms=12); ax.plot(r.C_only,yi,"o",color="C3" if r.p_only<0.05 else "C0",ms=7)
    ax.axvline(0.5,color="k",lw=0.6,ls=":"); ax.set_yticks(y); ax.set_yticklabels(df.index,fontsize=8); ax.set_title(title,fontsize=9)
    ax.set_xlabel("Cross-validated C-index (set alone)",fontsize=8); ax.set_xlim(0.38,0.78)
fig.suptitle("Published NSCLC miRNA sets vs random sets of the same size",fontsize=9)
plt.tight_layout(rect=[0,0,1,0.93]); plt.savefig(os.path.join(RES,"Figure1_published_vs_random_sets.png"),dpi=200)
