# SET-UP --------------------------------------------------------------
# Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",       # Misc. data manipulation and plotting
              "here",            # Managing file paths
              "clusterProfiler") # Enrichment testing
pacman::p_load(char = packages)

# Source script with functions
source("scripts/kegg_fun.R")

# Define the input files
# GFF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gff/data/>
DE_file <- here("results/DE/DE_all.txt")  # File with DE results
kegg_map_file <- file.path("results/kegg/kegg_map.txt")
kegg_descr_file <- file.path("results/kegg/kegg_pathways.txt")

# Define the output files
outdir <- here("results/kegg")
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile <- here(outdir, "kegg_agglomByTime.tsv")

# Settings - focal tissues and contrasts
tissues <- c("ab", "ht", "lrt")
contrasts <- c("lvp06_vir06", "lvp24_vir24")


# PREP FOR GO ANALYSIS ---------------------------------------------------------
# KEGG info
kegg_map <- read_tsv(kegg_map_file, show_col_types = FALSE)
kegg_descr <- read_tsv(kegg_descr_file, show_col_types = FALSE)

# Prep DE results
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


# RUN KEGG ANALYSIS ------------------------------------------------------------
combs <- expand_grid(times = c("t06", "t24"), DE_dirs = c("up", "down", "both"))
kegg_res_raw <- map2_dfr(.x = combs$times, .y = combs$DE_dirs, .f = run_kegg_time,
                         DE_res = DE_res, kegg_map = kegg_map)

kegg_res <- kegg_res_raw |> 
  left_join(kegg_descr, by = c("ID" = "pathway")) |>
  select(contrast,
         DE_direction = direction,
         padj = p.adjust,
         n_DE_in_cat = Count,
         category = ID,
         description,
         gene_ids = geneID) |>
  filter(padj < 0.05,
         n_DE_in_cat > 1,
         !is.na(description))

write_tsv(kegg_res, outfile)
