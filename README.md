# 16S_QIIME2
Snakemake based 16S QIIME2 pipeline 

## input
- Multiplexed R1/R2 read pairs
- QIIME2 compatible mapping files
  - Tab delimited
  - First two columns should be SampleID (or #SampleID) and BarcodeSequence

## demultiplex step
- Multiplexed -> demultiplexed fastq(.gz) files
- QIIME2 manifest file
- Total read count summary
