#!/usr/bin/env bash
# Reproduces the analysis. Put the raw files in data/raw first (see data/README.md).
# Settings used for the reported results: NREP=5, NRAND=300 (full run takes roughly 1 hour on a laptop).
set -e
cd "$(dirname "$0")/work"
export DATA_RAW=${DATA_RAW:-../data/raw} RESULTS_DIR=${RESULTS_DIR:-../results}
NREP=${NREP:-5}; NRAND=${NRAND:-300}; export NR=$NRAND
python ../scripts/01_prepare_tcga.py
for C in LUAD LUSC; do python ../scripts/03_audit_tcga.py $NREP $NRAND $C; done
for C in LUAD LUSC; do for K in 7 4; do python ../scripts/04_extra_null_sizes.py $C $K $NRAND; done; done
python ../scripts/05_rescore_all_signatures.py
python ../scripts/06_univariate_adjusted.py
python ../scripts/07_external_validation_GSE16025.py
python ../scripts/08_make_figure1.py
