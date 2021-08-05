# 16S_QIIME2
This is a Snakemake based 16S QIIME2 pipeline.

## Installation
To install, we assume you already have installed `Miniconda3 (4.7.10+)` (https://docs.conda.io/en/latest/miniconda.html)
- Clone the repository:
```bash
git clone https://github.com/junglee0713/16S_QIIME2.git
```
- Create a conda environment:
```bash
cd 16S_QIIME2
conda env create -f environment.yml
```
- The following software also need to be installed:
  - `dnabc` (https://github.com/PennChopMicrobiomeProgram/dnabc)
  - `unassigner` (https://github.com/kylebittinger/unassigner)

To run the pipeline, activate the envrionment (currently based on **QIIME2 2019.7**) by entering
`source activate qiime2-2019.7`
(Make sure to install `dnabc` and `unassigner` in the `qiime2-2019.7` envrionment)

## Required input files for the pipeline
To run the pipeline, we need
- Multiplexed R1/R2 read pairs (Undetermined_S0_L001_R1_001.fastq.gz, Undetermined_S0_L001_R2_001.fastq.gz), and
- QIIME2 compatible mapping file
  - Tab delimited
  - The first two columns should be `SampleID` (or `#SampleID`) and `BarcodeSequence`

## How to run
- Create a project directory, e.g. `/home/leej39/16S_QIIME2/test` and put the mapping file, e.g. `test_mapping_file.tsv` in the project directory
- Edit `config.yml` so that it suits your project. In particular,
  - **all: project**: path to the project directory, e.g. `/home/leej39/16S_QIIME2/test`
  - **all: mux_dir**: the direcotry containing multiplexed R1/R2 read pairs, e.g. `/home/leej39/16S_QIIME2/test/multiplexed_fastq` 
  - **all: mapping**: the name of mapping file, e.g. `test_mapping_file.tsv`
- To run the pipeline, activate the envrionment by entering `source activate qiime2-2019.7`, `cd` into `16S_QIIME2` and execute e.g.
```bash
snakemake \
    --configfile path/to/config.yml \
    --keep-going \
    --latency-wait 90 \
    --notemp \
    --printshellcmds
```
- When submitting jobs using `qsub`, you may run e.g. `bash run_snakemake.bash path/to/config_test.yml`
  - `bash dryrun_snakemake.bash path/to/config_test.yml` for dryrun
  - `bash unlock_snakemake.bash path/to/config_test.yml` for unlocking
  
## Intermediate steps and corresponding input/output

### Demultiplexing
#### Input
- Multiplexed R1/R2 read pairs
- QIIME2 compatible mapping file
#### Output
- Demultiplexed fastq(.gz) files
- Total read count summary (tsv)
- QIIME2 compatible manifest file (csv)

### QIIME2 import
#### Input
- QIIME2 compatible manifest file
- Demultiplexed fastq files
#### Output
- QIIME2 PairedEndSequencesWithQuality artifact and corresponding visualization
- QIIME2-generated demultiplexing stats

### DADA2 denoise
#### Input
- QIIME2 PairedEndSequencesWithQuality artifact
#### Output
- Feature table (QIIME2 artifact, tsv)
- Representative sequences (QIIME2 artifact, fasta)

### Taxonomy classification
#### Input
- Representative sequences 
#### Output
- Taxonomy classification table (QIIME2 artifact, tsv) 

### Tree building
#### Input
- Representative sequences 
#### Output
- Aligned sequence
- Masked (aligned) sequence
- Unrooted tree
- Rooted tree

### Diversity calculation
#### Input
- Rooted tree
#### Output
- Various QIIME2 diversity metric artifacts
- Faith phylogenetic diversity vector (tsv)
- Weighted/unweighted UniFrac distance matrices (tsv)

### Unassigner
#### Input
- Representative sequences (fasta)
#### Output
- Unassigner output (tsv) for species level classification of representative sequences

### dada2_species
#### Input
- Representative sequences (fasta)
#### Output
- Dada2 species assignments (tsv)
- Dada2 Raw data for loading in R (RData format)

### vsearch
#### Input
- Representative sequences (fasta)
#### Output
- Vsearch report (tsv) customized to be like BLAST results (see config.yml)
- Vsearch list of representative sequences that aligned (fasta)