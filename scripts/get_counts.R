## Load packages
if (!require(pacman)) install.packages("pacman")
pacs <- c("tidyverse", "here", "readxl", "janitor")
p_load(char = pacs)

## Define input files
abd_in <- here("results/trex/3084-Abd-DEgenes.xlsx")
ht_in <- here("results/trex/3084-HT-DEgenes.xlsx")
lrt_in <- here("results/trex/3084-LRT-DEgenes.xlsx")
lrt12_in <- here("results/trex/2048-Ferdinand-mosquito-VectorbaseV49ref-DEgenes.xlsx")

## Define output files
outdir <- here("results/counts")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile <- here(outdir, "counts_all.txt")

## Read and process file 1 (old data - sugar level 12%)
geneIDs <- read_xlsx(lrt12_in, sheet = "DE genes table", range = c("C5:C50000")) %>%
  rename(gene_id = GeneID)

lrt12_counts <- read_xlsx(lrt12_in, sheet = "DE genes table", range = c("DO5:EL50000")) %>%
  cbind(geneIDs, .) %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
  mutate(sample = sub("norm\\.", "", sample),
         # ID this tissue by a lowercase 'l', Al241 => Akh + LRT12% + 24h + rep1
         treatment = substr(sample, 1, 2),
         sample = sub("^([A-Z])", "\\1l", sample),
         tissue = "lrt12",
         sugar = "sug12",
         treatment = case_when(
           treatment == "A2" ~ "akh24",
           treatment == "A6" ~ "akh06",
           treatment == "L2" ~ "lvp24",
           treatment == "L6" ~ "lvp06",
           treatment == "V2" ~ "vir24",
           treatment == "V6" ~ "vir06",
         ))

## Function to read the other 3 files (new data - sugar level 3%)
process_xl <- function(xls_file, tissue) {
  geneIDs <- read_xlsx(xls_file, sheet = "DE genes table", range = c("C5:C50000")) %>%
    rename(gene_id = GeneID)
  
  read_xlsx(xls_file, sheet = "DE genes table", range = c("DO5:EF50000")) %>%
    cbind(geneIDs, .) %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "count") %>%
    mutate(sample = sub("norm\\.E2\\.", "", sample),
           tissue = tissue,
           sugar = "sug03",
           treatment = paste0(substr(sample, 1, 1), substr(sample, 3, 3)),
           treatment = case_when(
             treatment == "A2" ~ "akh24",
             treatment == "A6" ~ "akh06",
             treatment == "L2" ~ "lvp24",
             treatment == "L6" ~ "lvp06",
             treatment == "V2" ~ "vir24",
             treatment == "V6" ~ "vir06",
           ))
}

## Read the sugar level 3% files
ab_counts <- process_xl(abd_in, "ab")
ht_counts <- process_xl(ht_in, "ht")
lrt_counts <- process_xl(lrt_in, "lrt")

## Combine all dataframes
all_counts <- bind_rows(ab_counts, ht_counts, lrt_counts, lrt12_counts)

## Write output file
write_tsv(all_counts, outfile)
