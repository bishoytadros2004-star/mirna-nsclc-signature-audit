# Step 4: add random-set nulls for extra set sizes (we ran sizes 7 and 4). Usage: python 04_extra_null_sizes.py COHORT SIZE [NR]
import sys, json, numpy as np, pandas as pd, warnings
from sklearn.model_selection import KFold
from statsmodels.duration.hazard_regression import PHReg
from lifelines.utils import concordance_index
warnings.filterwarnings("ignore")
COH=sys.argv[1]; K=int(sys.argv[2]); NR=int(sys.argv[3]) if len(sys.argv)>3 else 300
rng=np.random.default_rng(77+K)
d=pd.read_pickle(f"{COH}.pkl").dropna(subset=['stage','age'])
mir=[c for c in d.columns if c.startswith('hsa-')]; idx={m:i for i,m in enumerate(mir)}
T=d['OS.time'].values/365.25; E=d['OS'].values.astype(int)
clin=d[['stage','age','male']].values.astype(float); X_all=d[mir].values
pool=[m for m in mir if d[m].median()>=2]
J=json.load(open(f"result_{COH}.json")); base=J["res"]["baseline_C"]
def cv(cols,bi,seed):
    Xm=X_all[:,[idx[c] for c in cols]]; X=np.hstack([clin,Xm]) if bi else Xm; cs=[]
    for tr,te in KFold(5,shuffle=True,random_state=seed).split(X):
        mu=X[tr].mean(0); sd=X[tr].std(0); keep=sd>1e-9
        Z=(X[:,keep]-mu[keep])/sd[keep]
        try:
            r=PHReg(T[tr],Z[tr],status=E[tr]).fit(disp=0); cs.append(concordance_index(T[te],-(Z[te]@r.params),E[te]))
        except Exception: pass
    return float(np.mean(cs))
o=[];p=[]
for j in range(NR):
    cols=list(rng.choice(pool,K,replace=False)); o.append(cv(cols,False,10+j)); p.append(cv(cols,True,500+j)-base)
J["null"][str(K)]={"only":o,"delta":p}
json.dump(J,open(f"result_{COH}.json","w")); print("done",COH,K)
