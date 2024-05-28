# GO wrapper function
run_GO_wrap <- function(fcontrast, direction, DE_res, GO_map, gene_lens) {
  DE_vec <- DE2vec(DE_res, fcontrast, direction)
  GO_df <- run_GO(fcontrast, direction, DE_vec, GO_map, gene_lens)
}

# GO wrapper function when aggregating DE results by time
run_GO_wrap_time <- function(ftime, direction, DE_res, GO_map, gene_lens) {
  DE_vec <- DE2vec_time(DE_res, ftime, direction)
  GO_df <- run_GO(fcontrast = ftime, direction, DE_vec, GO_map, gene_lens)
}

# Create named vector of DE genes
DE2vec <- function(DE_res, fcontrast, direction = "both") {
    if (direction == "up") DE_res <- DE_res %>% filter(lfc > 0)
  if (direction == "down") DE_res <- DE_res %>% filter(lfc < 0)
  
  DE_foc <- DE_res %>%
    filter(contrast == fcontrast,
           !is.na(padj)) %>%   # Exclude genes with NA adj-p-val - not tested
    mutate(sig = ifelse(padj < 0.05, 1, 0))

  contrast_vector <- DE_foc$sig
  names(contrast_vector) <- DE_foc$gene_id

  return(contrast_vector)
}

# Create named vector for DE results aggregated by timepoint
DE2vec_time <- function(DE_res, ftime, direction = "both") {
  if (direction == "up") DE_res <- DE_res %>% filter(lfc > 0)
  if (direction == "down") DE_res <- DE_res %>% filter(lfc < 0)
  
  DE_foc <- DE_res |> 
    filter(time == ftime, !is.na(padj)) |>
    arrange(gene_id, padj) |>
    slice_head(n = 1, by = gene_id) |>
    mutate(sig = ifelse(padj < 0.05, 1, 0))
  
  contrast_vector <- DE_foc$sig
  names(contrast_vector) <- DE_foc$gene_id
  
  return(contrast_vector)
}

# Function to run GO analysis
run_GO <- function(fcontrast, direction = "both",
                   DE_vec, GO_map, gene_lens) {

  if(sum(DE_vec > 0)) {
    # Remove rows from gene length df not in the DE_vec
    gene_lens_final <- gene_lens %>% filter(gene_id %in% names(DE_vec))

    # Remove elements from DE_vec not among the gene lengths
    DE_vec_final <- DE_vec[names(DE_vec) %in% gene_lens_final$gene_id]

    # Probability weighting function based on gene lengths
    pwf <- suppressMessages(
      nullp(DEgenes = DE_vec_final,
            bias.data = gene_lens_final$gene_length,
            plot.fit = FALSE)
    )

    # Run GO test
    GO_df <- suppressMessages(
      goseq(pwf = pwf, gene2cat = GO_map,
            method = "Wallenius", use_genes_without_cat = FALSE)
    )

    # Process GO results
    DEGs <- names(DE_vec)[DE_vec == 1]
    
    GO_df <- GO_df %>%
      filter(numDEInCat > 0) %>%    # P-adjustment only for genes that were actually tested
      mutate(padj = p.adjust(over_represented_pvalue, method = "BH"),
             contrast = fcontrast,
             direction = direction,
             numDE = length(DEGs)) %>%
      select(contrast, direction, padj, numDEInCat, numInCat, numDE,
             category, ontology, description = term)
    
    GO_df <- GO_map %>%
      filter(gene_id %in% DEGs,
             go_term %in% GO_df$category) %>%
      group_by(go_term) %>%
      summarize(DE_in_cat = paste(gene_id, collapse = ",")) %>%
      full_join(GO_df, ., by = c("category" = "go_term"))
    
    nsig <- GO_df |> filter(padj < 0.05, numDEInCat > 1) |>  nrow()
    cat("\n--------------\nRan comparison for:", fcontrast, direction, "\n")
    cat("# Number of significant DE genes:", sum(DE_vec), "\n")
    cat("# Number of significant GO categories:", nsig, "\n")

    return(GO_df)
    
  } else {
    
    cat("# No significant DE genes\n")
  
  }
}

# Prep df for GO plot
prep_goplot <- function(df, contrasts, tissues, directions = "both") {
  
  # Subset to focal contrasts and tissues
  df <- df %>%
    filter(contrast %in% contrasts,
           tissue %in% tissues,
           direction %in% directions)
  
  # When using >1 tissue, use longer contrast ID
  if (length(tissues) > 1) {
    df <- df %>% mutate(contrast = contrast_full)
    contrasts <- unique(df$contrast_full)
  }
  
  vars <- c("category", "ontology", "description",
            "tissue", "contrast", "direction", "padj")
  
  df %>%
    select(any_of(vars)) %>%
    # Pivot wider and then longer to include all terms in all contrasts
    pivot_wider(names_from = contrast, values_from = padj) %>%
    pivot_longer(cols = any_of(contrasts),
                 names_to = "contrast", values_to = "padj") %>%
    left_join(df %>% select(contrast, tissue, category, direction, numDEInCat, sig),
              by = c("contrast", "tissue", "category", "direction")) %>%
    # No labels if not significant
    mutate(numDEInCat = ifelse(padj >= 0.05, NA, numDEInCat)) %>%
    mutate(contrast = sub("padj_", "", contrast),
           padj = ifelse(sig == FALSE, NA, padj),
           padj_log = -log10(padj)) %>% 
    # Only take GO categories with at least one significant contrast
    filter(category %in% (filter(., sig == TRUE) %>% pull(category))) %>%
    # Only take contrast with at least one significant category
    filter(contrast %in% (filter(., sig == TRUE) %>% pull(contrast))) %>%
    arrange(padj_log) %>% 
    mutate(description = str_trunc(description, width = 45),
           description = ifelse(is.na(description), category, description),
           description = fct_inorder(description))
}

# Heatmap-style plot for significant GO categories
goplot <- function(df, x_var = "contrast", facet_var = NULL,
                   type = "GO", label_count = TRUE,
                   title = NULL, xlabs = NULL,
                   ylabsize = 9, numDE_size = 1.5) {
  p <- ggplot(df) +
    aes(x = .data[[x_var]], y = description, fill = padj_log) +
    geom_tile(stat = "identity", size = 0.25, color = "grey80") +
    scale_fill_viridis_c(option = "D", na.value = "grey95") +
    labs(fill = "-log10\n(adj. p)", title = title) +
    scale_y_discrete(position = "right") +
    theme_minimal() +
    theme(legend.position = "left",
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(size = ylabsize))

  if(label_count == TRUE) {
    p <- p + geom_label(aes(label = numDEInCat),
                        fill = "grey95", size = numDE_size)
  }
  
  if (type == "GO") {
    if (is.null(facet_var)) {
      # ggforce::facet_col will keep tile heights constant
      p <- p +
        facet_col(vars(ontology), scales = "free_y", space = "free") +
        theme(strip.text = element_text(size = 10, face = "bold"))
    } else {
      p <- p +
        facet_grid(rows = vars(ontology),
                   cols = vars(!!sym(facet_var)),
                   scales = "free_y",
                   space = "free_y",
                   switch = "y")
    }
  } else if (!is.null(facet_var)) {
    p <- p +
      #facet_col(vars(!!sym(facet_var)),
      #          scales = "free_y", space = "free") +
      facet_wrap(vars(!!sym(facet_var)), nrow = 1) +
      theme(strip.text = element_text(size = 10, face = "bold"))
  }
  
  if (!is.null(xlabs)) p <- p + scale_x_discrete(labels = xlabs)

  return(p)
}

# GO dotplot
godotplot <- function(df, type = "GO") {
  
  #if (type == "GO") group_by <- "ontology" else group_by <- NULL 
  
  p <- ggdotchart(df,
             x = "description", y = "padj_log",
             color = "padj_log",
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             rotate = TRUE,                                # Rotate vertically
             #group = group_by,                             # Order by groups
             dot.size = 5,                                 # Large dot size
             label = "numDEInCat",                         # Add nr DE genes as dot labels
             font.label = list(color = "white", size = 9, vjust = 0.5),
             ggtheme = theme_bw()) +                       # ggplot2 theme
    labs(y = "-log10(adj. p-value)", x = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_color_viridis_c(option = "D", na.value = "grey95") +
    theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
          plot.title = element_text(size = 15, face = "bold"),
          strip.text.y = element_text(angle = 270, face = "bold"),
          strip.placement = "outside",
          axis.title.x = element_text(margin = margin(t = 0.5, b = 0.5, unit = "cm")),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 8),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank())
  
  if (type == "GO") {
    p <- p + facet_grid(ontology~contrast, space = "free", scales = "free")
  } else if (type == "KEGG") {
    p <- p + facet_wrap(vars(contrast), scales = "free_x", nrow = 1)
  }
  
  print(p)
}
