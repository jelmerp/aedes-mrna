# SET-UP --------------------------------------------------------------
# Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "goseq",           # GO statistics
              "ape")             # To read GFF
pacman::p_load(char = packages)

# Source script with functions
source("scripts/GO_fun.R")

# Define the input files
# GFF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gff/data/>
DE_file <- here("results/DE/DE_all.txt")  # File with DE results
GO_map_file <- here("results/GO/GO_map.txt")
gff_file <- here("data/refdata/VectorBase-49_AaegyptiLVP_AGWG.gff")

# Define the output files
outdir <- here("results", "GO")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
GO_res_file <- here(outdir, "GO_agglomByTime.tsv")

# Settings - focal tissues and contrasts
tissues <- c("ab", "ht", "lrt") # "lrt12" - removed on Laura's request 2023-12-18
contrasts <- c("lvp06_vir06", "lvp24_vir24")


# PREP FOR GO ANALYSIS ---------------------------------------------------------
# Get GO mappings from Vectorbase GAF file
GO_map <- as.data.frame(read_tsv(GO_map_file, show_col_types = FALSE))

# Extract gene lengths from GFF file
gene_lens <- read.gff(gff_file) %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub("ID=(.*);.*", "\\1", attributes),
         gene_length = end - start) %>%
  select(gene_id, gene_length) %>%
  arrange(gene_id)

# Create named vector with significant and non-significant genes
DE_res <- read_tsv(DE_file, show_col_types = FALSE) |>
  mutate(tissue = gsub(".*\\d(\\w+)", "\\1", treat_a),
         tissue = ifelse(sugar == "sug12", "lrt12", tissue),
         trt_a = gsub("(.*\\d+)\\w+", "\\1", treat_a),
         trt_b = gsub("(.*\\d+)\\w+", "\\1", treat_b),
         time = sub(".*(06|24)", "t\\1", trt_a),
         contrast = paste0(trt_a, "_", trt_b),
         contrast_full = paste0(contrast, "_", tissue)) |>
  # Only select focal contrasts and tissues!
  filter(contrast %in% contrasts,
         tissue %in% tissues) |>
  select(gene_id, lfc, padj, trt_a, trt_b, tissue, time, contrast_full) |>
  arrange(contrast_full, gene_id)


# RUN GO ANALYSIS --------------------------------------------------------------
combs <- expand_grid(times = c("t06", "t24"), DE_dirs = c("up", "down", "both"))
GO_res <- map2_dfr(.x = combs$times, .y = combs$DE_dirs, .f = run_GO_wrap_time,
                   DE_res, GO_map, gene_lens)

GO_sig <- GO_res |>
  filter(padj < 0.05, numDEInCat > 1, !is.na(description)) |>
  rename(DE_direction = "direction")

write_tsv(GO_sig, GO_res_file)
