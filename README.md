# 16S_QIIME2
This is a Snakemake based 16S QIIME2 pipeline.

## Software requirement
We assume that the following software are installed:
- QIIME2 version 2019.4 (https://docs.qiime2.org/2019.4)
- dnabc (https://github.com/PennChopMicrobiomeProgram/dnabc)
- unassigner (https://github.com/kylebittinger/unassigner)

## input
- Multiplexed R1/R2 read pairs
- QIIME2 compatible mapping files
  - Tab delimited
  - First two columns should be `SampleID` (or `#SampleID`) and `BarcodeSequence`

## demultiplex step
- Multiplexed -> demultiplexed fastq(.gz) files
- QIIME2 manifest file
- Total read count summary
