# SET-UP -----------------------------------------------------------------------
## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "clusterProfiler") # Enrichment testing
pacman::p_load(char = packages)

## Source script with functions
source("scripts/kegg_fun.R")

## Input files
infile_DE <- here("results", "DE", "DE_all.txt")
kegg_dir <- here("results", "kegg")
infile_kegg_map <- file.path(kegg_dir, "kegg_map.txt")
infile_kegg_descr <- file.path(kegg_dir, "kegg_pathways.txt")

## Output file
outfile_kegg_res <- here(kegg_dir, "kegg_all.txt")

## Read input files
DE_res <- read_tsv(infile_DE, show_col_types = FALSE) %>%
  mutate(tissue = gsub(".*\\d(\\w+)", "\\1", treat_a),
         tissue = ifelse(sugar == "sug12", "lrt12", tissue),
         treat_a = gsub("(.*\\d+)\\w+", "\\1", treat_a),
         treat_b = gsub("(.*\\d+)\\w+", "\\1", treat_b),
         contrast = paste0(treat_a, "_", treat_b, "_", tissue)) %>% 
  arrange(contrast, gene_id)

kegg_map <- read_tsv(infile_kegg_map, show_col_types = FALSE)

kegg_descr <- read_tsv(infile_kegg_descr, show_col_types = FALSE)


# KEGG ENRICHMENT TEST ---------------------------------------------------------
## Get all pairwise contrasts
contrasts <- unique(DE_res$contrast)

## Run enrichment test for all pairwise contrasts
kegg_both <- map_dfr(.x = contrasts, .f = run_kegg, "both", DE_res, kegg_map)
kegg_up <- map_dfr(.x = contrasts, .f = run_kegg, "up", DE_res, kegg_map)
kegg_down <- map_dfr(.x = contrasts, .f = run_kegg, "down", DE_res, kegg_map)

kegg_res <- bind_rows(kegg_both, kegg_up, kegg_down) %>%
  left_join(kegg_descr, by = c("ID" = "pathway")) %>% 
  select(contrast,
         direction,
         padj = p.adjust,
         count = Count,
         category = ID,
         description,
         gene_ids = geneID)

## Write output file
write_tsv(kegg_res, outfile_kegg_res)
