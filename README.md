# 16S_QIIME2
This is a Snakemake based 16S QIIME2 pipeline.

## Software requirement
We assume that the following software are installed:
- `QIIME2 version 2019.4` (https://docs.qiime2.org/2019.4)
- `dnabc` (https://github.com/PennChopMicrobiomeProgram/dnabc)
- `unassigner` (https://github.com/kylebittinger/unassigner)

## Input
To run the pipeline, we need
- Multiplexed R1/R2 read pairs (Undetermined_S0_L001_R1_001.fastq.gz, Undetermined_S0_L001_R2_001.fastq.gz), and
- QIIME2 compatible mapping file
  - Tab delimited
  - The first two columns should be `SampleID` (or `#SampleID`) and `BarcodeSequence`

## Intermediate steps and correspondign input/output
### Demultiplexing
#### Input
- Multiplexed R1/R2 read pairs
- QIIME2 compatible mapping file
#### Output
- Demultiplexed fastq(.gz) files
- Total read count summary
- QIIME2 compatible manifest file

### QIIME2 import
#### Input
- QIIME2 compatible manifest file
- Demultiplexed fastq files
#### Output
- QIIME2 PairedEndSequencesWithQuality artifact and corresponding visualization
- QIIME2-generated demultiplexing stats
