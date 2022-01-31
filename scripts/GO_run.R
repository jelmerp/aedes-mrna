# SET-UP --------------------------------------------------------------
## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "goseq",           # GO statistics
              "ape")             # To read GFF
pacman::p_load(char = packages)

## Source script with functions
source("scripts/GO_fun.R")

## Input files
### GAF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gaf/>
### GFF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gff/data/>
DE_file <- here("results/DE/DE_all.txt")  # File with DE results
gaf_file <- here("refdata", "VectorBase-49_AaegyptiLVP_AGWG_GO.gaf")
gff_file <- here("refdata", "VectorBase-49_AaegyptiLVP_AGWG.gff")

## Output files
outdir <- here("results", "GO")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
GO_res_file <- here(outdir, "GO_all.txt")


# PREP FOR GO ANALYSIS ---------------------------------------------------------
## Create named vector with significant and non-significant genes
DE_res <- read_tsv(DE_file, show_col_types = FALSE) %>%
  mutate(tissue = gsub(".*\\d(\\w+)", "\\1", treat_a),
         tissue = ifelse(sugar == "sug10", "lrt10", tissue),
         treat_a = gsub("(.*\\d+)\\w+", "\\1", treat_a),
         treat_b = gsub("(.*\\d+)\\w+", "\\1", treat_b),
         contrast = paste0(treat_a, "_", treat_b, "_", tissue)) %>% 
  arrange(contrast, gene_id)

## Get GO mappings from Vectorbase GAF file
GO_map <- read_tsv(gaf_file, skip = 1,
                   col_names = FALSE, show_col_types = FALSE) %>%
  select(gene_id = X2, go_term = X5) %>%
  arrange(gene_id, go_term) %>%
  data.frame(.)  # Needs to be a df and not a tibble for GOseq

## Extract gene lengths from GFF file
gene_lens <- read.gff(gff_file) %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub("ID=(.*);.*", "\\1", attributes),
         gene_length = end - start) %>%
  select(gene_id, gene_length) %>%
  arrange(gene_id)


# RUN GO ANALYSIS --------------------------------------------------------------
GO_res <- map_dfr(.x = unique(DE_res$contrast), .f = run_GO_wrap,
                  DE_res, GO_map, gene_lens)
write_tsv(GO_res, GO_res_file)
