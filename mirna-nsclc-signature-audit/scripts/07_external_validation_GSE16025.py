# Step 7: independent squamous cohort GSE16025 (61 tumours). Needs GSE16025 series matrix + GPL5106 platform file in data/raw.
import os, io, gzip, re, json, numpy as np, pandas as pd, warnings
from importlib import import_module
from sklearn.model_selection import KFold
from statsmodels.duration.hazard_regression import PHReg
from lifelines.utils import concordance_index
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
RAW=os.environ.get("DATA_RAW","../data/raw")
def pick(*names):
    for n in names:
        p=os.path.join(RAW,n)
        if os.path.exists(p): return p
    raise FileNotFoundError("missing one of: "+", ".join(names))
# --- series matrix
meta={}; 
with gzip.open(pick("GSE16025_series_matrix.txt.gz","GSE16025_series_matrix_txt.gz"),"rt") as f:
    lines=f.read().split("\n")
gsm=src=None
for l in lines:
    if l.startswith("!Sample_geo_accession"): gsm=[x.strip('"') for x in l.split("\t")[1:]]
    if l.startswith("!Sample_source_name_ch1"): src=[x.strip('"') for x in l.split("\t")[1:]]
    if l.startswith("!Sample_characteristics_ch1"):
        v=[x.strip('"') for x in l.split("\t")[1:]]; k=v[0].split(":")[0].strip(); meta[k]=[(x.split(":",1)[1].strip() if ":" in x else "") for x in v]
b=lines.index("!series_matrix_table_begin"); e=lines.index("!series_matrix_table_end")
tab=pd.read_csv(pd.io.common.StringIO("\n".join(lines[b+1:e])),sep="\t",index_col=0)
tab.index=[str(i).strip('"') for i in tab.index]; tab.columns=[c.strip('"') for c in tab.columns]
clin=pd.DataFrame(meta,index=gsm); clin['src']=src
clin=clin[clin.src.str.contains("SCC")]
clin['T']=pd.to_numeric(clin['survival_time'])/12.0; clin['E']=pd.to_numeric(clin['death']).astype(int)
clin['age']=pd.to_numeric(clin['age']); clin['stage']=clin['stage'].str.upper().str.extract(r'(I+V?)')[0].map({'I':1,'II':2,'III':3,'IV':4})
X=tab[clin.index].T
print("tumours",len(clin),"events",int(clin.E.sum()),"expression range",float(np.nanmin(X.values)),float(np.nanmax(X.values)),"NaN cells",int(np.isnan(X.values).sum()))
# --- probe -> miRNA base name
_ln=[];_in=False
with gzip.open(pick("GPL5106_family.soft.gz","GPL5106_family_soft.gz"),"rt",errors="ignore") as _f:
    for _l in _f:
        if _l.startswith("!platform_table_begin"): _in=True; continue
        if _l.startswith("!platform_table_end"): break
        if _in: _ln.append(_l)
P=pd.read_csv(io.StringIO("".join(_ln)),sep="\t",low_memory=False).set_index("ID")
def base(n):
    n=str(n).lower().replace("hsa-","")
    n=re.sub(r'(\*|-5p|-3p)$','',n)
    return n
names={pid:base(P.loc[pid,'Human Mature Name']) for pid in X.columns if pid in P.index and isinstance(P.loc[pid,'Human Mature Name'],str)}
F=X[list(names)].copy(); F.columns=[names[c] for c in F.columns]
F=F.T.groupby(level=0).mean().T
F=F.loc[:,[c for c in F.columns if c.startswith(("mir-","let-"))]]
F=F.fillna(F.mean())
print("miRNA features after collapsing probes:",F.shape[1])
SIGB={"Li2017":["mir-21","mir-30d"],"Gallach2017":["mir-21","mir-188"],"MacDonagh2021":["mir-30a","mir-30b","mir-30c","mir-34a","mir-4286"],
 "Hu2010":["mir-486","mir-30d","mir-1","mir-499"],"Guo2016":["mir-1","mir-30d","mir-221","mir-486"],"Zhang2019":["mir-191","mir-28","mir-145","mir-328","mir-18a"],
 "Liu2017":["mir-23b","mir-10b","mir-21"],"Zhou2020":["mir-942","mir-601"],"Xue2020":["mir-1228","mir-181a"],"Yang2019":["mir-146b","mir-29c"],
 "Raponi2009":["mir-146b"],"Yu2008":["let-7a","mir-221","mir-137","mir-372","mir-182"],"Effusion2012":["mir-93","mir-100","mir-134","mir-151","mir-345"]}
T=clin['T'].values; E=clin['E'].values; cl=clin[['stage','age']].values.astype(float)
NREP=int(os.environ.get('NREP_EXT','10'))
def cv(cols,bi,seed,reps=NREP):
    M=F[cols].values if cols else np.zeros((len(F),0)); Xc=np.hstack([cl,M]) if bi else M; cs=[]
    for r in range(reps):
        for tr,te in KFold(5,shuffle=True,random_state=seed*100+r).split(Xc):
            mu=Xc[tr].mean(0); sd=Xc[tr].std(0); keep=sd>1e-9
            if keep.sum()==0: continue
            Z=(Xc[:,keep]-mu[keep])/sd[keep]
            try:
                res=PHReg(T[tr],Z[tr],status=E[tr]).fit(disp=0); cs.append(concordance_index(T[te],-(Z[te]@res.params),E[te]))
            except Exception: pass
    return float(np.mean(cs))
rng=np.random.default_rng(16025)
pool=list(F.columns[F.mean()>=F.mean().quantile(0.5)])
base_C=cv([],True,1); print("baseline (stage+age) CV C-index %.3f ; expressed pool %d"%(base_C,len(pool)))
rows=[]; NR=int(os.environ.get('NR','300')); nulls={}
for n,mirs in SIGB.items():
    have=[m for m in mirs if m in F.columns]; k=len(have)
    miss=[m for m in mirs if m not in F.columns]
    if k==0: rows.append(dict(set=n,k=0,missing=";".join(miss))); continue
    if k not in nulls:
        o=[];p=[]
        for j in range(NR):
            cols=list(rng.choice(pool,k,replace=False)); o.append(cv(cols,False,10+j,2)); p.append(cv(cols,True,500+j,2)-base_C)
        nulls[k]=(np.array(o),np.array(p))
    co=cv(have,False,2); cp=cv(have,True,3)-base_C; no,nd=nulls[k]
    rows.append(dict(set=n,k=k,missing=";".join(miss),C_only=co,null_med=float(np.median(no)),p_only=float((np.sum(no>=co)+1)/(NR+1)),
                     delta=cp,null_delta_med=float(np.median(nd)),p_delta=float((np.sum(nd>=cp)+1)/(NR+1))))
df=pd.DataFrame(rows)
ok=df.k>0
df.loc[ok,'q_only']=multipletests(df[ok].p_only,method='fdr_bh')[1]; df.loc[ok,'q_delta']=multipletests(df[ok].p_delta,method='fdr_bh')[1]
df.to_csv(os.path.join(os.environ.get("RESULTS_DIR","../results"),"external_validation_GSE16025_LUSC.csv"),index=False)
pd.set_option('display.width',200)
print(df.round(3).to_string(index=False))
print("nominal p_only<0.05:",int((df.p_only<0.05).sum()),"of",int(ok.sum()),"; FDR<0.05:",int((df.q_only<0.05).sum()))
print("nominal p_delta<0.05:",int((df.p_delta<0.05).sum()),"; positive delta:",int((df.delta>0).sum()),"; max delta %.3f"%df.delta.max())
