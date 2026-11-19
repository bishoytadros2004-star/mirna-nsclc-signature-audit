# Step 1: build tumour-only miRNA + survival + clinical tables. Run from the work/ folder.
import os, pandas as pd, numpy as np
RAW=os.environ.get("DATA_RAW","../data/raw")
def pick(*names):
    for n in names:
        p=os.path.join(RAW,n)
        if os.path.exists(p): return p
    raise FileNotFoundError("missing one of: "+", ".join(names))
out={}
for coh in ["LUAD","LUSC"]:
    m=pd.read_csv(pick(f"TCGA-{coh}.mirna.tsv.gz",f"TCGA-{coh}_mirna_tsv.gz"),sep="\t",index_col=0)
    s=pd.read_csv(pick(f"TCGA-{coh}.survival.tsv.gz",f"TCGA-{coh}_survival_tsv.gz"),sep="\t")
    c=pd.read_csv(pick(f"TCGA-{coh}.GDC_phenotype.tsv.gz",f"TCGA-{coh}_clinical_tsv.gz"),sep="\t",low_memory=False)
    # primary tumour samples only (-01), one per patient
    tum=[x for x in m.columns if x.split('-')[3].startswith('01')]
    mt=m[tum].T; mt['patient']=[x[:12] for x in mt.index]
    mt=mt[~mt.patient.duplicated()].set_index('patient')
    s['patient']=s['_PATIENT']; s=s.drop_duplicates('patient').set_index('patient')
    c['patient']=c['sample'].str[:12]; c=c.drop_duplicates('patient').set_index('patient')
    stage=c['ajcc_pathologic_stage.diagnoses'].astype(str).str.upper()
    def st(x):
        if 'IV' in x: return 4
        if 'III' in x: return 3
        if 'II' in x: return 2
        if 'I' in x: return 1
        return np.nan
    clin=pd.DataFrame({'stage':stage.map(st),'age':pd.to_numeric(c['age_at_index.demographic'],errors='coerce'),
                       'male':(c['gender.demographic']=='male').astype(float)})
    d=s[['OS.time','OS']].join(clin,how='inner').join(mt.drop(columns=[],errors='ignore'),how='inner')
    d=d[(d['OS.time']>0)&d['OS'].notna()]
    print(coh,"tumours with miRNA:",mt.shape[0],"with survival:",len(d),"events:",int(d.OS.sum()),"median follow-up days:",d['OS.time'].median(),
          "stage NaN:",int(d.stage.isna().sum()),"age NaN:",int(d.age.isna().sum()))
    d.to_pickle(f"{coh}.pkl")
    mir=mt.columns
    ex=(mt.median()>=2)
    print("  precursors total",len(mir),"median log2RPM>=2:",int(ex.sum()),"mean>=1:",int((mt.mean()>=1).sum()))
