# SET-UP -----------------------------------------------------------------------
# Install/load packages
dyn.load("/fs/ess/PAS0471/jelmer/software/GLPK/lib/libglpk.so.40", local=FALSE)
library(tidyverse)
library(here)
library(clusterProfiler)

# Input files
infile_DE <- here("results/DE/DE_all.txt")
kegg_dir <- here("results/kegg")
infile_kegg_map <- here(kegg_dir, "KEGG_map.txt")
infile_kegg_descr <- here(kegg_dir, "KEGG_pathways.txt")
# Lists of genes DE in both our and Alfonso-Parra's study (emailed by Laura on 2024-05-22)
DE_6hpm_file <- "results/DE/DE_overlap_Alfonso-Parra_6hpm.txt"
DE_24hpm_file <- "results/DE/DE_overlap_Alfonso-Parra_24hpm.txt"

# Output file
outfile_kegg_res <- here(kegg_dir, "kegg_studyoverlap_20240529.tsv")

# Read input files
DE_6hpm <- readLines(DE_6hpm_file)
DE_24hpm <- readLines(DE_24hpm_file)
kegg_map <- read_tsv(infile_kegg_map, show_col_types = FALSE)
kegg_descr <- read_tsv(infile_kegg_descr, show_col_types = FALSE)
allgenes <- read_tsv(infile_DE, show_col_types = FALSE) |>
  distinct(gene_id) 

# KEGG ENRICHMENT TEST ---------------------------------------------------------
res_6hpm <- as.data.frame(enricher(DE_6hpm, TERM2GENE = kegg_map)) |>
  mutate(contrast = "6hpm")
res_24hpm <- as.data.frame(enricher(DE_24hpm, TERM2GENE = kegg_map)) |>
  mutate(contrast = "24hpm")

kegg_res <- bind_rows(res_6hpm, res_24hpm) |>
  left_join(kegg_descr, by = c("ID" = "pathway")) |>
  select(contrast,
         padj = p.adjust,
         count = Count,
         category = ID,
         description,
         gene_ids = geneID)

# Write output file
write_tsv(kegg_res, outfile_kegg_res)
