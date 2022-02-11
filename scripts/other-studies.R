# SETUP ------------------------------------------------------------------------
## Install/load packages
if(! "pacman" %in% installed.packages()) install.packages("pacman")
packages <- c("tidyverse",        # Misc. data manipulation and plotting
              "here",             # Managing file paths
              "readxl")           # Reading Excel files
pacman::p_load(char = packages)

## Define input files
indir_comp <- here("refdata/other-studies")

alfpar_basename <- "Alfonso-Parra2016_LRT_updated-from-Amaro.xlsx"
alonso_basename <- "41598_2019_52268_Alonso et al 2019 HT_10%.xlsx"
pascini_basename <- "12864_2020_6543_MOESM5_ESM. Pascini et al., 2020xlsx.xlsx"
camargo_basename <- "Camargo_S3_41598_2020_71904_MOESM3_ESM.xlsx"
amaro_basename <- "Amaro BMCGenomics 2021 Ae MAG_Saline 24h DEG Summary.xlsx"
  
alonso_file <- here(indir_comp, alonso_basename)
alfpar_file <- here(indir_comp, alfpar_basename)
pascini_file <- here(indir_comp, pascini_basename)
camargo_file <- here(indir_comp, camargo_basename)
amaro_file <- here(indir_comp, amaro_basename)

## Define output file
outdir <- "results/other-studies"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
all_out <- here(outdir, "other-studies_combined.txt")


# READ XLS FILES ---------------------------------------------------------------
## Read and process Alonso et al 2019 results
alonso_raw <- read_excel(alonso_file, sheet = "results_VRG_INS",
                         skip = 1, col_types = c("text", rep("numeric", 6)))
colnames(alonso_raw)[1] <- "gene_id"
alonso <- alonso_raw %>%
  mutate(gene_id = sub("(.*)-.*", "\\1", gene_id), study = "alonso") %>%
  rename(lfc = log2FoldChange) %>%
  select(gene_id, lfc, padj, study) %>%
  group_by(gene_id) %>%
  arrange(padj) %>%
  slice(1) # ONLY TAKE THE TRANSCRIPT WITH THE LOWEST P-ADJ

## Read and process Alfonso-Parra et al 2016 results
alfpar_24h <- read_excel(alfpar_file, sheet = "DE_significant_genes 24 hpm") %>% 
  select(gene_id = GeneID, lfc = logFC, padj = FDR) %>%
  mutate(study = "alfpar", hpm = 24)
alfpar_6h <- read_excel(alfpar_file, sheet = "DE_significant_genes 6 hpm") %>% 
  select(gene_id = GeneID, lfc = logFC, padj = FDR) %>%
  mutate(study = "alfpar", hpm = 6)

## Read and process Pascini et al 2020 results
pascini <- read_excel(pascini_file, sheet = "Aedes-CDS") %>% 
  select(gene_id = "Best match to AEGY5.2 database...31",
         #gene_id = "Best match to AEGY-CDS database",
         rpkm_virg = "RPKM Virgin...119",
         rpkm_ins = "RPKM Ins...120",
         padj = "FDR...341") %>%
  mutate(lfc = log2(rpkm_virg / rpkm_ins),
         study = "pascini") %>%
  select(gene_id, lfc, padj, study) %>%
  filter(!is.na(gene_id), gene_id != "A") %>%
  group_by(gene_id) %>%
  arrange(padj) %>%
  slice(1) # ONLY TAKE THE TRANSCRIPT WITH THE LOWEST P-ADJ

## Read and process Camargo at al 2020 results
camargo <- read_excel(camargo_file, sheet = "DE transcrips") %>%
  filter(Trt %in% c("NBF6", "NBF24")) %>%                         # Non-blood-fed at 6 and 24hpm
  mutate(hpm = as.integer(sub("NBF", "", Trt))) %>% 
  select(gene_id = VB_ID, lfc = logFC, padj = FDR, hpm) %>%
  mutate(study = "camargo") %>%
  filter(gene_id != "NA")
camargo_6h <- camargo %>% filter(hpm == 6)
camargo_24h <- camargo %>% filter(hpm == 24)

## Read and process Amaro at al 2021 results
amaro1 <- read_excel(amaro_file, sheet = "Ae_MAG vs Ae_none DEG")
amaro_control <- read_excel(amaro_file, sheet = "Ae_Saline vs Ae_none DEG")
amaro <-  anti_join(amaro1, amaro_control, by = "GeneID") %>%
  select(gene_id = GeneID, lfc = logFC, padj = FDR, hpm = time) %>%
  mutate(hpm = as.integer(sub("^Hr", "", hpm)),
         study = "amaro")

# COMBINE DFs ------------------------------------------------------------------
all <- bind_rows(alfpar_6h, alfpar_24h, alonso, amaro, camargo, pascini)
write_tsv(all, all_out)
