# Step 3: set-level cross-validated C-index + random-set null. Usage: python 03_audit_tcga.py NREP NRAND COHORT
import sys, time, json, warnings, numpy as np, pandas as pd
from sklearn.model_selection import KFold
from statsmodels.duration.hazard_regression import PHReg
from lifelines.utils import concordance_index
from importlib import import_module
SIGS=import_module('02_signatures').SIGS
warnings.filterwarnings("ignore")
NREP=int(sys.argv[1]); NRAND=int(sys.argv[2]); COH=sys.argv[3]
rng=np.random.default_rng(2026)
d=pd.read_pickle(f"{COH}.pkl").dropna(subset=['stage','age'])
mir=[c for c in d.columns if c.startswith('hsa-')]
T=d['OS.time'].values/365.25; E=d['OS'].values.astype(int)
clin=d[['stage','age','male']].values.astype(float)
pool=[m for m in mir if d[m].median()>=2]
X_all=d[mir].values; idx={m:i for i,m in enumerate(mir)}
def cv(cols, base_included, reps=NREP, seed=0):
    Xm=X_all[:,[idx[c] for c in cols]] if cols else np.zeros((len(d),0))
    X=np.hstack([clin,Xm]) if base_included else Xm
    cs=[]
    for r in range(reps):
        kf=KFold(5,shuffle=True,random_state=seed*100+r)
        for tr,te in kf.split(X):
            mu=X[tr].mean(0); sd=X[tr].std(0); sd[sd==0]=1
            Z=(X-mu)/sd
            try:
                res=PHReg(T[tr],Z[tr],status=E[tr]).fit(disp=0)
                lp=Z[te]@res.params
                cs.append(concordance_index(T[te],-lp,E[te]))
            except Exception:
                pass
    return float(np.mean(cs)) if cs else np.nan
res={"cohort":COH,"n":len(d),"events":int(E.sum())}
t0=time.time()
base=cv([],True,seed=1)   # clinical only (stage, age, sex)
res["baseline_C"]=base
res["sets"]={}
sizes=set()
for n,src,typ,mirs in SIGS:
    have=[m for m in mirs if m in idx]
    sizes.add(len(have))
    res["sets"][n]={"src":src,"type":typ,"k":len(have),"only_C":cv(have,False,seed=2),"plus_C":cv(have,True,seed=3)}
    res["sets"][n]["delta"]=res["sets"][n]["plus_C"]-base
print(COH,"real sets done in %.0fs"%(time.time()-t0),flush=True)
null={}
for k in sorted(sizes):
    o=[];p=[]
    for j in range(NRAND):
        cols=list(rng.choice(pool,k,replace=False))
        o.append(cv(cols,False,reps=1,seed=10+j)); p.append(cv(cols,True,reps=1,seed=500+j)-base)
    null[k]={"only":o,"delta":p}
    print(COH,"null size",k,"done %.0fs"%(time.time()-t0),flush=True)
for n in res["sets"]:
    k=res["sets"][n]["k"]; s=res["sets"][n]
    s["p_only"]=float((np.sum(np.array(null[k]["only"])>=s["only_C"])+1)/(NRAND+1))
    s["p_delta"]=float((np.sum(np.array(null[k]["delta"])>=s["delta"])+1)/(NRAND+1))
    s["null_only_median"]=float(np.median(null[k]["only"])); s["null_delta_median"]=float(np.median(null[k]["delta"]))
json.dump({"res":res,"null":{str(k):v for k,v in null.items()}},open(f"result_{COH}.json","w"))
print("FINISHED",COH,"%.0fs"%(time.time()-t0),flush=True)
