#!/bin/bash
snakemake -j 80 \
	--configfile config.yml \
	--cluster-config cluster.json \
	-w 90 \
	--notemp \
	-p \
	-c \
	"qsub -cwd -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"
