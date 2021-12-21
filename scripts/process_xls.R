## Load packages
if (!require(pacman)) install.packages("pacman")
pacs <- c("tidyverse", "here", "readxl", "janitor")
p_load(char = pacs)

## Define input files
abd_in <- here("trex-docs/3084-Abd-DEgenes.xlsx")
ht_in <- here("trex-docs/3084-HT-DEgenes.xlsx")
lrt_in <- here("trex-docs/3084-LRT-DEgenes.xlsx")

## Define output files
outdir <- here("results/DE")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
outfile_all <- here(outdir, "DE_all.txt")
outfile_sig <- here(outdir, "DE_sig.txt")
outfile_idx <- here(outdir, "gene_info.txt")

## Function to create a df for DE results for 1 pairwise comparison
get1comp <- function(treat_a, treat_b,
                     xls, sheet = "DE genes table", range,
                     idx) {
  de <- read_xlsx(xls, sheet = sheet, range = range)%>%
    bind_cols(idx, .) %>%
    drop_na() %>%
    clean_names() %>%
    mutate(treat_a = treat_a,
           treat_b = treat_b,
           contrast = paste0(treat_a, "_", treat_b),
           padj = as.numeric(padj)) %>%
    select(gene_id, treat_a, treat_b, lfc = log2_fc, pvalue, padj, contrast)
}

## Function to extract DE results for all comparisons in 1 Excel file
get1file <- function(xls, tissue) {
  idx <- read_xlsx(xls, sheet = "DE genes table", range = c("C5:H50000"))
  
  a6 <- paste0("akh06", tissue)
  l6 <- paste0("lvp06", tissue)
  v6 <- paste0("vir06", tissue)
  a24 <- paste0("akh24", tissue)
  l24 <- paste0("lvp24", tissue)
  v24 <- paste0("vir24", tissue)
    
  av6 <- get1comp(a6, v6, xls, range = "M5:W50000", idx = idx)
  lv6 <- get1comp(l6, v6, xls, range = "X5:AH50000", idx = idx)
  la6 <- get1comp(l6, a6, xls, range = "AI5:AS50000", idx = idx)
  
  av24 <- get1comp(a24, v24, xls, range = "AT5:BD50000", idx = idx)
  lv24 <- get1comp(l24, v24, xls, range = "BE5:BO50000", idx = idx)
  la24 <- get1comp(l24, a24, xls, range = "BP5:BZ50000", idx = idx)
  
  time_a <- get1comp(a24, a6, xls, range = "CA5:CK50000", idx = idx)
  time_l <- get1comp(l24, l6, xls, range = "CL5:CV50000", idx = idx)
  time_v <- get1comp(v24, v6, xls, range = "CW5:DG50000", idx = idx)
  
  all <- bind_rows(av6, lv6, la6,
                   av24, lv24, la24,
                   time_a, time_l, time_v)
  return(all)
}

## Create df's with DE genes
abd <- get1file(abd_in, "ab")
ht <- get1file(ht_in, "ht")
lrt <- get1file(lrt_in, "lrt")

all <- bind_rows(abd, ht, lrt)
sig <- filter(all, padj < 0.05)

## Create df with gene info
idx <- read_xlsx(abd_in, sheet = "DE genes table", range = c("C5:H50000")) %>%
  rename(gene_id = GeneID)

## Write output files
write_tsv(all, outfile_all)
write_tsv(sig, outfile_sig)
write_tsv(idx, outfile_idx)
