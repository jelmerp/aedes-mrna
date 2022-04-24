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
GO_map_file <- here(outdir, "GO_map.txt")

# PREP FOR GO ANALYSIS ---------------------------------------------------------
## Create named vector with significant and non-significant genes
DE_res <- read_tsv(DE_file, show_col_types = FALSE) %>%
  mutate(tissue = gsub(".*\\d(\\w+)", "\\1", treat_a),
         tissue = ifelse(sugar == "sug12", "lrt12", tissue),
         treat_a = gsub("(.*\\d+)\\w+", "\\1", treat_a),
         treat_b = gsub("(.*\\d+)\\w+", "\\1", treat_b),
         contrast = paste0(treat_a, "_", treat_b, "_", tissue)) %>% 
  arrange(contrast, gene_id)

## SET CONTRASTS - don't use AKH and don't use across-time comps
CONTRASTS <- unique(DE_res$contrast)[!grepl("akh|06.*24|24.*06", unique(DE_res$contrast))]

## Get GO mappings from Vectorbase GAF file
GO_map <- read_tsv(gaf_file, skip = 1,
                   col_names = FALSE, show_col_types = FALSE) %>%
  select(gene_id = X2, go_term = X5) %>%
  arrange(gene_id, go_term) %>%
  distinct() %>% 
  data.frame(.)  # Needs to be a df and not a tibble for GOseq
write_tsv(GO_map, GO_map_file)

## Extract gene lengths from GFF file
gene_lens <- read.gff(gff_file) %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub("ID=(.*);.*", "\\1", attributes),
         gene_length = end - start) %>%
  select(gene_id, gene_length) %>%
  arrange(gene_id)


# RUN GO ANALYSIS --------------------------------------------------------------
GO_both <- map_dfr(CONTRASTS, run_GO_wrap, "both", DE_res, GO_map, gene_lens)
GO_up <- map_dfr(CONTRASTS, run_GO_wrap, "up", DE_res, GO_map, gene_lens)
GO_down <- map_dfr(CONTRASTS, run_GO_wrap, "down", DE_res, GO_map, gene_lens)
GO_res <- bind_rows(GO_both, GO_up, GO_down)
write_tsv(GO_res, GO_res_file)
