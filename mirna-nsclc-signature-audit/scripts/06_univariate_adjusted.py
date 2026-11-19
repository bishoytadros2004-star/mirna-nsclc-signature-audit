# Step 6: each miRNA from the published sets, Cox model adjusted for stage, age, sex (TCGA). Writes results/univariate_adj.csv
import os, numpy as np, pandas as pd, warnings
from importlib import import_module
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
SIGS=import_module("02_signatures").SIGS
mirs=sorted({m for _,_,_,l in SIGS for m in l}); rows=[]
for coh in ["LUAD","LUSC"]:
    d=pd.read_pickle(f"{coh}.pkl").dropna(subset=["stage","age"])
    for m in mirs:
        if m not in d.columns or d[m].std()<1e-6: continue
        x=d[[m,"OS.time","OS","stage","age","male"]].copy(); x[m]=(x[m]-x[m].mean())/x[m].std()
        try:
            r=CoxPHFitter().fit(x,"OS.time","OS").summary.loc[m]
            rows.append(dict(cohort=coh,miRNA=m,HR=float(np.exp(r["coef"])),lo=float(np.exp(r["coef lower 95%"])),hi=float(np.exp(r["coef upper 95%"])),p=float(r["p"])))
        except Exception as e: print("fail",coh,m,str(e)[:60])
df=pd.DataFrame(rows)
df["q"]=df.groupby("cohort").p.transform(lambda p: multipletests(p,method="fdr_bh")[1])
out=os.environ.get("RESULTS_DIR","../results"); os.makedirs(out,exist_ok=True)
df.to_csv(os.path.join(out,"univariate_adj.csv"),index=False); print(df.groupby("cohort").p.apply(lambda p:(p<0.05).sum()))
