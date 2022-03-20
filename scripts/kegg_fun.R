run_kegg <- function(fcontrast, DE_res, kegg_map) {
  
  DE_genes <- DE_res %>%
    filter(padj < 0.05, contrast == fcontrast) %>%
    pull(gene_id)
  
  cat("## Contrast:", fcontrast, " // Nr DE genes: ", length(DE_genes))
  if (length(DE_genes) >  1) {
    kegg_res <- enricher(DE_genes, TERM2GENE = kegg_map) %>%
      as.data.frame(.)
    cat(" // Nr enriched pathways:", nrow(kegg_res), "\n")
  
    if (nrow(kegg_res) > 0) {
      kegg_res$contrast <- fcontrast
      row.names(kegg_res) <- NULL
      return(kegg_res)
    }
  } else {
    cat("\n")
  }
}

## Get the NCBI gene ID
geneID_lookup <- function(geneID) {
  res <- entrez_search(db = "gene", term = geneID)$ids
  print(res)
  return(res)
}

## Get the genes belonging to a certain KEGG pathway
get_pw_genes <- function(pathway_id) {
  print(pathway_id)
  
  pw <- keggGet(pathway_id)
  if (is.null(pw[[1]]$GENE)) return(NA)
  
  pw2 <- pw[[1]]$GENE[c(TRUE, FALSE)]
  pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
  
  return(pw2)
}
