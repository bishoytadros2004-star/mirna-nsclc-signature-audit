# Step 5: final scoring of all 13 sets against the saved random-set nulls; writes results/audit_results_by_signature.csv
import os, json, numpy as np, pandas as pd, warnings
from sklearn.model_selection import KFold
from statsmodels.duration.hazard_regression import PHReg
from lifelines.utils import concordance_index
from statsmodels.stats.multitest import multipletests
from importlib import import_module
SIGS=import_module('02_signatures').SIGS
warnings.filterwarnings("ignore")
NREP=5; rows=[]
for COH in ["LUAD","LUSC"]:
    d=pd.read_pickle(f"{COH}.pkl").dropna(subset=['stage','age'])
    mir=[c for c in d.columns if c.startswith('hsa-')]; idx={m:i for i,m in enumerate(mir)}
    T=d['OS.time'].values/365.25; E=d['OS'].values.astype(int)
    clin=d[['stage','age','male']].values.astype(float); X_all=d[mir].values
    J=json.load(open(f"result_{COH}.json")); null=J["null"]; base=J["res"]["baseline_C"]
    def cv(cols,base_inc,seed):
        Xm=X_all[:,[idx[c] for c in cols]] if cols else np.zeros((len(d),0))
        X=np.hstack([clin,Xm]) if base_inc else Xm
        cs=[]; fails=0
        for r in range(NREP):
            for tr,te in KFold(5,shuffle=True,random_state=seed*100+r).split(X):
                mu=X[tr].mean(0); sd=X[tr].std(0); keep=sd>1e-9
                if keep.sum()==0: fails+=1; continue
                Z=(X[:,keep]-mu[keep])/sd[keep]
                try:
                    res=PHReg(T[tr],Z[tr],status=E[tr]).fit(disp=0)
                    cs.append(concordance_index(T[te],-(Z[te]@res.params),E[te]))
                except Exception: fails+=1
        return (float(np.mean(cs)) if cs else np.nan), fails
    for n,src,typ,mirs in SIGS:
        have=[m for m in mirs if m in idx]
        zero=[m for m in have if d[m].std()<1e-9]
        co,f1=cv(have,False,2); cp,f2=cv(have,True,3)
        k=len(have); no=np.array(null[str(k)]["only"]); nd=np.array(null[str(k)]["delta"])
        rows.append(dict(cohort=COH,set=n,type=typ,k=k,zero_var=len(zero),failed_folds=f1+f2,C_only=co,null_med=float(np.nanmedian(no)),
                         p_only=float((np.sum(no>=co)+1)/(len(no)+1)),delta=cp-base,null_delta_med=float(np.nanmedian(nd)),
                         p_delta=float((np.sum(nd>=cp-base)+1)/(len(nd)+1))))
df=pd.DataFrame(rows)
assert df.C_only.notna().all()
df['q_only']=df.groupby('cohort').p_only.transform(lambda p: multipletests(p,method='fdr_bh')[1])
df['q_delta']=df.groupby('cohort').p_delta.transform(lambda p: multipletests(p,method='fdr_bh')[1])
df.to_csv(os.path.join(os.environ.get("RESULTS_DIR","../results"),"audit_results_by_signature.csv"),index=False)
pd.set_option('display.width',220)
print(df[['cohort','set','k','zero_var','failed_folds','C_only','null_med','p_only','q_only','delta','p_delta','q_delta']].round(3).to_string(index=False))
print("\nnominal p_only<0.05: %d of %d; FDR q<0.05: %d"%((df.p_only<0.05).sum(),len(df),(df.q_only<0.05).sum()))
print("nominal p_delta<0.05: %d of %d; FDR q<0.05: %d"%((df.p_delta<0.05).sum(),len(df),(df.q_delta<0.05).sum()))
print("positive incremental C: %d of %d; max delta %.3f"%((df.delta>0).sum(),len(df),df.delta.max()))
