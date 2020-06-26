#!/usr/bin/env Rscript

library(tidyverse)
library(Biostrings)
library(dada2)

#from config.yml

rdp_species <- snakemake@config[["dada2"]][["species_training_set"]]

#from dada2.rules

my_seqs <- readDNAStringSet(snakemake@input[[1]])

#main

seqs_df <- tibble(OTU_ID = names(my_seqs),
                  sequence = as.character(my_seqs))

genus.species <- assignSpecies(seqs_df$sequence, rdp_species)

genspec_df <- genus.species %>%
    as.data.frame() %>%
    rownames_to_column(var = "sequence") %>%
    inner_join(seqs_df, by = "sequence") %>%
    select(OTU_ID, Genus, Species)

#from dada2.rules and targets.rules

write_tsv(x = genspec_df, path = snakemake@output)
