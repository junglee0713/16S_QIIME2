#!/usr/bin/env bash

#SBATCH --mem=2G
#SBATCH -n 1
#SBATCH --export=ALL
#SBATCH --mail-user=tanesc@chop.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --no-requeue
#SBATCH -t 12:00:00
#SBATCH --output=slurm_%x_%j.out

#Uncomment the next two lines if you want to 'qsub' this script
#source ~/.bashrc #needed to make "conda" command to work
#conda activate qiime2-snakemake

set -xeuo pipefail

if [ $# -ne 1 ]; then
    echo "Usage: bash $0 PATH_TO_CONFIG"
    exit 1
fi

CONFIG_FP=$1

snakemake \
    --jobs 100 \
    --configfile ${CONFIG_FP} \
    --cluster-config cluster.json \
    --keep-going \
    --latency-wait 90 \
    --notemp \
    --printshellcmds \
    --cluster \
    "sbatch --mem=10G -t 12:00:00 -n 4 --export=ALL --no-requeue"
