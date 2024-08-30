# 2024-08-30, Jelmer Poelstra, R version 4.4.0

# Load packages
library(tidyverse)

# Define in- and output files
infile <- "results/counts/counts_all.txt"
outfile <- "results/counts/count_matrix.tsv"

# Read the input files
count_df <- read_tsv(infile, show_col_types = FALSE)

# Change the long-format count dataframe into a wide count matrix-like format
count_mat <- count_df |> 
  separate_wider_position(treatment, widths = c(mating = 3, hpm = 2)) |> 
  filter(mating != "akh") |>
  mutate(
    hpm = paste0("hpm", hpm),
    tissue = sub("lrt12", "lrt", tissue),
    mating = ifelse(mating == "vir", "unmated", "mated"),
    sample_full = paste(sample, tissue, mating, sugar, hpm, sep = "_")
    ) |>
  select(gene_id, sample_full, count) |>
  filter(!is.na(gene_id)) |> 
  pivot_wider(id_cols = "gene_id", names_from = "sample_full", values_from = "count")

# Write the output files
write_tsv(count_mat, outfile)
