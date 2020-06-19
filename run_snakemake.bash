#!/bin/bash
set -xeuo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: bash $0 PATH_TO_CONFIG"
    exit 1
fi

CONFIG_FP=$1

source ~/.bashrc
conda activate qiime2-snakemake

snakemake \
    --jobs 100 \
    --configfile ${CONFIG_FP} \
    --cluster-config cluster.json \
    --keep-going \
    --latency-wait 90 \
    --notemp \
    --printshellcmds \
    --cluster \
    "qsub -cwd -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"
