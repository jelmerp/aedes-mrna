# 2024-05-28 -- Run GO analysis on genes overlapping between our study and the
#               Alfonso-Parra study.
# Running this in R 4.3.0

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

# Input files
# GAF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gaf/>
# GFF file from: <https://vectorbase.org/vectorbase/app/downloads/release-49/AaegyptiLVP_AGWG/gff/data/>
DE_file <- here("results/DE/DE_all.txt")  # File with DE results, used to get background list of genes
gaf_file <- here("data/ref/aedes/VectorBase-49_AaegyptiLVP_AGWG_GO.gaf")
gff_file <- here("data/ref/aedes/VectorBase-49_AaegyptiLVP_AGWG.gff")

# Lists of genes DE in both our and Alfonso-Parra's study (emailed by Laura on 2024-05-22)
DE_6hpm_file <- "results/DE/DE_overlap_Alfonso-Parra_6hpm.txt"
DE_24hpm_file <- "results/DE/DE_overlap_Alfonso-Parra_24hpm.txt"

# Output files
outdir <- here("results/GO")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
GO_res_file <- here(outdir, "GO_studyoverlap_20240529.tsv")


# PREP FOR GO ANALYSIS ---------------------------------------------------------
# Read files with DE genes
DE_6hpm <- readLines(DE_6hpm_file)
DE_24hpm <- readLines(DE_24hpm_file)

# Get vector of background universe of genes
allgenes <- read_tsv(DE_file, show_col_types = FALSE) |>
  distinct(gene_id) |>
  mutate(sig = 0)

# Create named vectors with significant and non-significant genes
DE_df_6hpm <- allgenes |> mutate(sig = ifelse(gene_id %in% DE_6hpm, 1, 0))
DE_vec_6hpm <- DE_df_6hpm$sig
names(DE_vec_6hpm) <- DE_df_6hpm$gene_id
#table(DE_vec_6hpm)
#length(DE_6hpm)

DE_df_24hpm <- allgenes |> mutate(sig = ifelse(gene_id %in% DE_24hpm, 1, 0))
DE_vec_24hpm <- DE_df_24hpm$sig
names(DE_vec_24hpm) <- DE_df_24hpm$gene_id
#table(DE_vec_24hpm)
#length(DE_24hpm)

# Get GO mappings from Vectorbase GAF file
GO_map <- read_tsv(gaf_file, skip = 1,
                   col_names = FALSE, show_col_types = FALSE) %>%
  select(gene_id = X2, go_term = X5) %>%
  arrange(gene_id, go_term) %>%
  distinct() %>% 
  data.frame(.)  # Needs to be a df and not a tibble for GOseq

# Extract gene lengths from GFF file
gene_lens <- read.gff(gff_file) %>%
  filter(type == "gene") %>%
  mutate(gene_id = sub("ID=(.*);.*", "\\1", attributes),
         gene_length = end - start) %>%
  select(gene_id, gene_length) %>%
  arrange(gene_id)


# RUN GO ANALYSIS --------------------------------------------------------------
GO_res_6hpm <- run_GO(fcontrast = "6hpm", direction = "NA",
                      DE_vec = DE_vec_6hpm, GO_map = GO_map, gene_lens = gene_lens) |>
  select(!c(direction, numInCat, numDE))
GO_res_24hpm <- run_GO(fcontrast = "24hpm", direction = "NA",
                      DE_vec = DE_vec_24hpm, GO_map = GO_map, gene_lens = gene_lens) |>
  select(!c(direction, numInCat, numDE))

GO_res <- bind_rows(GO_res_6hpm, GO_res_24hpm) |>
  filter(padj < 0.05)
write_tsv(GO_res, GO_res_file)
