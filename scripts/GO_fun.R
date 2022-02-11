run_GO_wrap <- function(fcontrast, DE_res, GO_map, gene_lens) {
  DE_vec <- get_DE_vec(fcontrast, DE_res)
  GO_df <- run_GO(fcontrast, DE_vec, GO_map, gene_lens)
}

## Create named vector of DE genes
get_DE_vec <- function(fcontrast, DE_res) {
  DE_foc <- DE_res %>%
    filter(contrast == fcontrast,
           !is.na(padj)) %>%   # Exclude genes with NA adj-p-val - not tested
    mutate(sig = ifelse(padj < 0.05, 1, 0))

  contrast_vector <- DE_foc$sig
  names(contrast_vector) <- DE_foc$gene_id

  return(contrast_vector)
}

## Function to run GO analysis
run_GO <- function(fcontrast, DE_vec, GO_map, gene_lens) {

  if(sum(DE_vec > 0)) {
    ## Remove rows from gene length df not in the DE_vec
    gene_lens_final <- gene_lens %>% filter(gene_id %in% names(DE_vec))

    ## Remove elements from DE_vec not among the gene lengths
    DE_vec_final <- DE_vec[names(DE_vec) %in% gene_lens_final$gene_id]

    ## Probability weighting function based on gene lengths
    pwf <- nullp(
      DEgenes = DE_vec_final,
      bias.data = gene_lens_final$gene_length,
      plot.fit = FALSE
    )

    ## Run GO test
    GO_df <- goseq(pwf = pwf, gene2cat = GO_map, method = "Wallenius")

    ## Process GO results
    GO_df <- GO_df %>%
      filter(numDEInCat > 0) %>%    # P-adjustment only for genes that were actually tested
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             contrast = fcontrast) %>%
      select(contrast, padj, numDEInCat, numInCat, category, ontology,
             description = term)

    n_sig <- nrow(filter(GO_df, padj < 0.05))
    cat("\n--------------\nRan comparison for:", fcontrast, "\n")
    cat("## Number of significant DE genes:", sum(DE_vec), "\n")
    cat("## Number of significant GO categories:", n_sig, "\n")

    return(GO_df)
  } else {
    cat("## No significant DE genes\n")
  }
}

## Heatmap-style plot for significant GO categories
goplot <- function(GO_df, x_var = "contrast",
                   title = NULL, xlabs = NULL, ylabsize = 9) {
  p <- GO_df %>%
    mutate(description = str_trunc(description, width = 45),
           description = ifelse(is.na(description), category, description)) %>%
    ggplot(aes_string(x_var, "description", fill = "padj")) +
    geom_tile(stat = "identity", size = 0.25, color = "grey40") +
    scale_fill_distiller(palette = "Reds", na.value = "grey90") +
    labs(fill = "adjusted\np-value",
         title = title) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = ylabsize))

  if(!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)

  return(p)
}

## Prep df for GO plot
prep_goplot <- function(GO_df, contrasts, tissues) {
  
  missing <- contrasts[!contrasts %in% GO_df$contrast]
  if (length(missing) > 0) GO_df <- GO_df %>% add_row(contrast = missing)
  
  GO_df <- GO_df %>% filter(contrast %in% contrasts, tissue %in% tissues)
  if (length(tissues) > 1) GO_df <- GO_df %>% mutate(contrast = contrast_full)
  
  GO_df %>%
    select(contrast, padj, category, description) %>%
    pivot_wider(names_from = contrast, values_from = padj) %>%
    pivot_longer(cols = -c(category, description),
                 names_to = "contrast", values_to = "padj") %>%
    mutate(padj = ifelse(padj >= 0.05, NA, padj)) %>% 
    filter(category %in% (filter(., padj < 0.05) %>% pull(category)))
}
