# GENES DE IN MULTIPLE STUDIES -------------------------------------------------
get_mult <- function(df) {
  mult_sig <- df %>%
    rowwise(gene_id, description, hpm) %>%
    mutate(`nr. signif.` = sum(c_across(contains("padj")) < 0.05, na.rm = TRUE)) %>%
    filter(padj_ours < 0.05, `nr. signif.` > 1) %>% 
    relocate(`nr. signif.`, .after = hpm) %>%
    arrange(-`nr. signif.`) %>%
    ungroup()
  
  ## Make the table print-friendly
  mult_sig <- mult_sig %>% rename(`gene ID` = gene_id) 
  colnames(mult_sig) <- sub("lfc_", "LFC: ", colnames(mult_sig))
  colnames(mult_sig) <- sub("padj_", "p: ", colnames(mult_sig))
  
  return(mult_sig)
}

## Create datatable
multsig_dt <- function(df) {
  
  df <- get_mult(df) %>% select(-hpm)
  
  p_cols <- colnames(df)[grep("p: ", colnames(df))] 
  
  df %>%
    #mutate(across(contains("p:"), f_sci)) %>%    # This doesn't work in combination with conditional formatting
    make_dt(caption = multsig_cap, pageLength = 10) %>%
    sigcolfmt(colnames = p_cols)
}

## Format significant columns with DT
sigcolfmt <- function(dt, colnames) {
  dt %>% formatStyle(
    columns = c(colnames),
    backgroundColor = styleInterval(0.05, c("lightgreen", "white"))
  )
}


# LFC PLOTS --------------------------------------------------------------------
prep_lfc <- function(df, study1, study2,
                     DE_levels = c("Neither", "x only", "y only", "Both")) {
  
  padj_study1 <- paste0("padj_", study1)
  padj_study2 <- paste0("padj_", study2)
  
  df %>%
    filter(!is.na(.data[[padj_study1]]),
           !is.na(.data[[padj_study2]])) %>%
    select(gene_id, description, contains(study1), contains(study2)) %>% 
    mutate(DE_in = case_when(
      .data[[padj_study1]] < 0.05 & .data[[padj_study2]] < 0.05 ~ "Both",
      .data[[padj_study1]] < 0.05 & .data[[padj_study2]] >= 0.05 ~ "x only",
      .data[[padj_study1]] >= 0.05 & .data[[padj_study2]] < 0.05 ~ "y only",
      .data[[padj_study1]] >= 0.05 & .data[[padj_study2]] >= 0.05 ~ "Neither"
    ),
    DE_in = droplevels(factor(DE_in, levels = DE_levels))
    )
}

plot_lfc_ <- function(
  df, study1, study2, interactive = FALSE,
  DE_levels = c("Neither", "x only", "y only", "Both"),
  DE_cols = c("grey60", colorblindr::palette_OkabeIto[c(1, 4, 3)]),
  studies_long = c("This study", "Alfonso-Parra", "Alonso", "Amaro",
                   "Camargo", "Pascini")
  ) {
  
  lfc_study1 <- paste0("lfc_", study1)
  lfc_study2 <- paste0("lfc_", study2)
  
  name_study1 <- studies_long[which(study_levels == study1)]
  name_study2 <- studies_long[which(study_levels == study2)]
  
  fDE_levels <- levels(df$DE_in)
  fcols <- DE_cols[match(fDE_levels, DE_levels)]
  
  flabs <- sub("^x", name_study1, fDE_levels)
  flabs <- sub("^y", name_study2, flabs)
  
  p <- df %>%
    ggplot(aes(x = .data[[lfc_study1]],
               y = .data[[lfc_study2]],
               color = DE_in,
               text = glue("Gene: {gene_id}
                           Description: {description}"))) +
    geom_point() +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_color_manual(values = fcols, labels = flabs) +
    labs(x = paste("LFC:", name_study1),
         y = paste("LFC:", name_study2),
         color = "DE in:") +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          aspect.ratio = 1)
  
  if (interactive == TRUE) p <- ggplotly(p, tooltip = "text")
  return(p)
}

plot_lfc <- function(df, study1, study2, interactive = FALSE) {
  
  message("## Study 1: ", study1, " / study 2: ", study2)
  
  p <- prep_lfc(df, study1, study2) %>%
    plot_lfc_(study1, study2, interactive)
  
  if (interactive == TRUE) return(p) else print(p)
}