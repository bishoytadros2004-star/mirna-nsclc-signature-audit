# Replication audit of published microRNA prognostic signatures in NSCLC

Code and result tables for the manuscript *"Do published microRNA prognostic signatures for non-small cell lung cancer outperform random miRNA sets? A replication audit in TCGA and an independent cohort"* https://github.com/bishoytadros2004-star/mirna-nsclc-signature-audit.git
## What this repository does
For 13 published NSCLC miRNA sets (`data/signatures.csv`) it asks whether each set carries more prognostic information than random miRNA sets of the same size, and whether it adds information to stage, age and sex, in TCGA-LUAD, TCGA-LUSC and the independent squamous cohort GSE16025. Per-set Cox models are refitted on the set's miRNAs (set-level test; published weights were not available for most sets) and scored by repeated 5-fold cross-validated C-index; the random-set null uses 300 random sets per set size.

## Main result (see manuscript for details and limitations)
Published sets performed like random sets in TCGA and added no meaningful information beyond clinical variables. Two sets reached FDR < 0.05 in GSE16025; both depend on miR-146b, whose association was reported from that cohort's original publication and which is not prognostic in TCGA.

## Layout
```
data/signatures.csv      the 13 sets, sources and precursor names
data/README.md           where to download the raw data (not included)
scripts/                 numbered analysis steps (01 to 08)
run_all.sh               runs everything in order
results/                 result tables and Figure 1 from the reported run
work/                    intermediate files (created when running)
```

## How to run
```bash
pip install -r requirements.txt
# download the raw files into data/raw (see data/README.md)
./run_all.sh              # full settings (about 1 hour)
NREP=1 NRAND=3 NREP_EXT=2 ./run_all.sh   # quick smoke test; p-values are then meaningless
```
Random seeds are fixed in the scripts. The tables in `results/` come from the original run of these scripts; the complete pipeline has been smoke-tested with small settings but a full end-to-end re-run to confirm that every number is reproduced has not been done yet [update this line after you re-run].

## Limitations
13 sets only (5 tissue-derived, 8 blood or effusion); precursor-level TCGA data; set-level test rather than published formulas; one small independent cohort that is the discovery cohort of one audited set; short TCGA follow-up. See the manuscript.

## Use of generative AI
[Insert the same statement as in the manuscript, for example: AI assistance was used for literature-search support, statistical analysis code and manuscript drafting; the authors verified all analyses, numbers and references.]

## Citation and license
Cite the manuscript. Code: MIT license. Result tables and figure: CC BY 4.0.
