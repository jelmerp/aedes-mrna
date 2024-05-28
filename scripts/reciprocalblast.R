# SETUP ------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(GO.db)   # BiocManager::install("GO.db")

# Define the input files
# 1. DIAMOND results with Dmel as db
todmel_file <- "results/ortho/diamond/to_aedes/diamond_out.tsv"
# 2. DIAMOND results with Aedes as db
toaedes_file <- "results/ortho/diamond/to_dmel/diamond_out.tsv"
# 3. GO annotation file from FlyBase (see run/prep.sh)
GO_file <- "data/ref/dmel/gene_association.fb"
# 4. Transcript-to-gene lookup file for Dmel (created by run/run_ortho.sh)
#    We need this because DIAMOND was run on proteins with transcript IDs,
#    while the GO annotations are linked to genes
dmel_trans2gene_file <- "results/ortho/annot/dmel_translens.tsv"

# Define the output files
outfile <- "results/ortho/RBBH_withGO.tsv"


# READ THE INPUTS --------------------------------------------------------------
to_dmel <- read_tsv(todmel_file, col_names = FALSE, comment = "#", show_col_types = FALSE) |>
  select(aedes_q = X1, dmel_s = X2)
to_aedes <- read_tsv(toaedes_file, col_names = FALSE, comment = "#", show_col_types = FALSE) |>
  select(dmel_q = X1, aedes_s = X2)
GO_raw <- read_tsv(GO_file, comment = "!", col_names = FALSE, show_col_types = FALSE)
trans2gene <- read_tsv(dmel_trans2gene_file, show_col_types = FALSE,
                       col_names = c("dmel_protein", "dmel_gene", "length")) |>
  select(dmel_protein, dmel_gene)


# GET RBBHs --------------------------------------------------------------------
# Get the reciprocal best blast hits (RBBH)
brh <- inner_join(to_dmel, to_aedes, by = join_by("dmel_s" == "dmel_q")) |>
  filter(aedes_q == aedes_s) |>
  # Aedes gene ID is simply the protein ID without the suffix:
  mutate(aedes_gene = sub("-.*", "", aedes_q)) |>
  select(aedes_gene, dmel_protein = dmel_s)

# Check that each gene is only present once in RBBH dataframe
stopifnot(!any(duplicated(brh$aedes_gene)))
stopifnot(!any(duplicated(brh$dmel_protein)))


# ADD GO ANNOTATIONS AND FINALIZE ----------------------------------------------
# Create a df with GOterm-to-description rows
go_info <- AnnotationDbi::select(GO.db,
                                 columns = c("GOID", "TERM", "ONTOLOGY"),
                                 keys = keys(GO.db, keytype = "GOID"),
                                 keytype = "GOID") |>
  dplyr::rename(GO_term = GOID, GO_description = TERM, GO_ontology = "ONTOLOGY")

# Process the GO table downloaded from FlyBase
#> For column info, see https://wiki.flybase.org/wiki/FlyBase:Gene_Ontology_(GO)_Annotation
#> "The evidence code for the GO annotation; one of IMP, IGI, IPI, ISS, IDA, IEP, IEA, TAS, NAS, ND, IC, RCA, HDA, HMP, HGI, HEP, IBA"
GO <- GO_raw |>
  dplyr::select(dmel_gene = X2, GO_qualifier = X4, GO_term = X5, GO_evidence = X7) |>
  left_join(go_info, by = "GO_term")

# Add the GO info to the BRH table
final_df <- left_join(brh, trans2gene, by = "dmel_protein") |>
  left_join(GO, by = "dmel_gene") |>
  dplyr::select(-dmel_protein)

# Write the output to file
write_tsv(final_df, outfile)
