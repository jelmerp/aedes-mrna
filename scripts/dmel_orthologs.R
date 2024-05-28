# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(here)

# Define input files
DE_file <- here("results/DE/DE_all.txt")
gene_file <- here("results/DE/gene_info.txt")
GO_annot_file <- here("results/orthologs/RBBH_withGO.tsv")

# Define output files
outfile <- here("results/orthologs/DE_withDmelGO.tsv")

# Settings - focal tissues and contrasts
tissues <- c("ab", "ht", "lrt", "lrt12")
contrasts <- c("lvp06_vir06", "lvp24_vir24")


# MAIN -------------------------------------------------------------------------
# Read and process input files
GO_annot <- read_tsv(GO_annot_file, show_col_types = FALSE)

gene_df <- read_tsv(gene_file, show_col_types = FALSE) |>
  select(gene_id, gene_description = description, gene_type = type) |>
  left_join(GO_annot, by = join_by("gene_id" == "aedes_gene"))

DE_res <- read_tsv(DE_file, show_col_types = FALSE) |>
  mutate(tissue = gsub(".*\\d(\\w+)", "\\1", treat_a),
         tissue = ifelse(sugar == "sug12", "lrt12", tissue),
         trt_a = gsub("(.*\\d+)\\w+", "\\1", treat_a),
         trt_b = gsub("(.*\\d+)\\w+", "\\1", treat_b),
         contrast = paste0(trt_a, "_", trt_b),
         contrast_full = paste0(contrast, "_", tissue)) |>
  # Only select focal contrasts and tissues!
  filter(contrast %in% contrasts,
         tissue %in% tissues,
         padj < 0.05) |>
  select(gene_id, lfc, padj, trt_a, trt_b, tissue, contrast_full)

# Merge the gene_df, DE_res, and GO_annot df's
final_df <- left_join(DE_res, gene_df, by = "gene_id",
                      relationship = "many-to-many") |>
  arrange(contrast_full, padj)
  
# Write output
write_tsv(final_df, outfile)
