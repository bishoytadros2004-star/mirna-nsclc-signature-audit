# Data (not stored in this repository)

The raw files are public and are not redistributed here. Download them into `data/raw/` (keep the file names).

## TCGA (UCSC Xena, GDC hub)
Go to https://xenabrowser.net/datapages/ and open the cohorts **GDC TCGA Lung Adenocarcinoma (LUAD)** and
**GDC TCGA Lung Squamous Cell Carcinoma (LUSC)**. For each, download:

| File (Xena name) | What it is |
|---|---|
| `TCGA-LUAD.mirna.tsv.gz` / `TCGA-LUSC.mirna.tsv.gz` | miRNA expression, precursor level, log2(RPM+1) |
| `TCGA-LUAD.survival.tsv.gz` / `TCGA-LUSC.survival.tsv.gz` | overall survival (OS, OS.time in days) |
| `TCGA-LUAD.GDC_phenotype.tsv.gz` / `TCGA-LUSC.GDC_phenotype.tsv.gz` | clinical phenotype (stage, age, sex); the analysis used the file called "clinical" or "phenotype" on the cohort page |

The scripts also accept the names with underscores (for example `TCGA-LUAD_mirna_tsv.gz`) that some browsers produce.

## GSE16025 (Gene Expression Omnibus)
- Series matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16025/matrix/GSE16025_series_matrix.txt.gz
- Platform GPL5106: https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL5nnn/GPL5106/soft/GPL5106_family.soft.gz

Download date and database versions used for the reported results: [add dates before release].
