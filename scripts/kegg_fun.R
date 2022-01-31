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
