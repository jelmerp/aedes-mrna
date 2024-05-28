run_kegg_time <- function(ftime,
                          direction = "both",
                          DE_res,
                          kegg_map) {
  
  if (direction == "up") DE_res <- DE_res %>% filter(lfc > 0)
  if (direction == "down") DE_res <- DE_res %>% filter(lfc < 0)
  
  DE_genes <- DE_res |> 
    filter(time == ftime) |> 
    arrange(gene_id, padj) |>
    slice_head(n = 1, by = gene_id) |>
    filter(padj < 0.05) |>
    pull(gene_id)
  
  cat("## Time:", ftime,
      " // Direction:", direction,
      " // Nr DE genes: ", length(DE_genes))
  
  if (length(DE_genes) >  1) {
    kegg_res <- enricher(DE_genes, TERM2GENE = kegg_map) %>%
      as.data.frame(.)
    cat(" // Nr enriched pathways:", nrow(kegg_res), "\n")
    
    if (nrow(kegg_res) > 0) {
      kegg_res$contrast <- ftime
      kegg_res$direction <- direction
      row.names(kegg_res) <- NULL
      return(kegg_res)
    }
  } else {
    cat("\n")
  }
}

run_kegg <- function(fcontrast, direction = "both", DE_res, kegg_map) {
  
  if (direction == "up") DE_res <- DE_res %>% filter(lfc > 0)
  if (direction == "down") DE_res <- DE_res %>% filter(lfc < 0)
  
  DE_genes <- DE_res %>%
    filter(padj < 0.05, contrast == fcontrast) %>%
    pull(gene_id)
  
  cat("## Contrast:", fcontrast,
      " // Direction:", direction,
      " // Nr DE genes: ", length(DE_genes))
  
  if (length(DE_genes) >  1) {
    kegg_res <- enricher(DE_genes, TERM2GENE = kegg_map) %>%
      as.data.frame(.)
    cat(" // Nr enriched pathways:", nrow(kegg_res), "\n")
  
    if (nrow(kegg_res) > 0) {
      kegg_res$contrast <- fcontrast
      kegg_res$direction <- direction
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
