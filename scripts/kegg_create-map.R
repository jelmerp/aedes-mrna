# TO DO: RESOLVE AAEL-IDs with 0 and with >1 matches to NCBI IDs?

# SET-UP -----------------------------------------------------------------------
## Source script with functions
source("scripts/kegg_fun.R")

## Load packages
if(!"pacman" %in% installed.packages()) install.packages("pacman")
library(pacman)
packages <- c("KEGGREST", "rentrez", "tidyverse", "here")
p_load(char = packages, install = TRUE)

## Define input files
genes_file <- here("results/DE/2048-Ferdinand-mosquito-VectorbaseV49ref-DEgenes.txt")

## Define output files
outdir <- "results/kegg"
kegg_map_file <- here(outdir, "kegg_map.txt")
pathway_df_file <- here(outdir, "kegg_pathways.txt")

lookup_raw_file <- here(outdir, "geneID_lookup_raw.rds")
lookup_file <- here(outdir, "geneID_lookup.rds")

## Get gene IDs
geneIDs <- read_delim(DE_file, delim = "\t") %>%
  arrange(GeneID) %>%
  pull(GeneID)


# GENE ID CONVERSION -----------------------------------------------------------
## Using the "rentrez" package
## link <- entrez_link(dbfrom="gene", id="5579271", db="all")
## AAEL015651

lookup_raw <- data.frame()
for(i in 1:length(geneIDs)) {
  geneID <- geneIDs[i]
  cat(i, geneID, "\n")
  geneID_NCBI <- entrez_search(db = "gene", term = geneID)$ids
  if (is_empty(geneID_NCBI)) geneID_NCBI <- NA
  lookup_raw <- rbind(lookup_raw, cbind(geneID, geneID_NCBI))
}
saveRDS(lookup_raw, lookup_raw_file)

lookup <- lookup_raw %>%
  drop_na() %>%
  group_by(geneID) %>%
  filter(n() == 1) %>% # Remove AAEL_IDs with multiple NCBI matches
  ungroup()
saveRDS(lookup, lookup_file)


# Get list of genes by KEGG pathway --------------------------------------------
## Check if Aedes is available - yes:
## org <- keggList("organism")
## org[grep("aedes", org[, 3], ignore.case = TRUE), ] # organism = "aag"

# Get KEGG pathways for Aedes:
pathway_list <- keggList("pathway", "aag")
pathway_codes <- sub("path:", "", names(pathway_list))
pathway_df <- data.frame(pathway = pathway_codes,
                         description = unname(pathway_list))
write_tsv(pathway_df, pathway_df_file)

# Pull all genes for each pathway
genes_by_pathway_raw <- sapply(pathway_codes, get_pw_genes)
genes_by_pathway <- tibble(pathway = names(genes_by_pathway_raw),
                           genes = genes_by_pathway_raw) %>%
  unnest(cols = genes)
colnames(genes_by_pathway) <- c("KEGG_pathway", "geneID_NCBI")


# Create KEGG map --------------------------------------------------------------
kegg_map <- merge(genes_by_pathway, lookup, by = "geneID_NCBI") %>%
  select(kegg_pathway, geneID) %>%
  arrange(kegg_pathway)

write_tsv(kegg_map, kegg_map_file)
