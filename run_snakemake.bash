#!/bin/bash

snakemake -j 80 \
	--configfile /scr1/users/leej39/16S_QIIME2/config.yml \
	--cluster-config /scr1/users/leej39/16S_QIIME2/cluster.json \
	-w 90 \
	--notemp \
	-p \
	-c \
	"qsub -cwd -r n -V -l h_vmem={cluster.h_vmem} -l mem_free={cluster.mem_free} -pe smp {threads}"
