## Function to make a kable table
make_kable <- function(df, caption = NULL) {
  df %>%
    kable(caption = caption,
          format.args = list(big.mark = ",")) %>%
    kable_styling(full_width = FALSE, bootstrap_options = "striped")
}

## Function to make an exportable datatable
make_dt <- function(df, numr_cols = "auto", caption = NULL,
                    filter = filter, pageLength = 10,
                    simple_mode = FALSE) {
  
  if (simple_mode == TRUE) {
    dom <- "t"
    paging <- FALSE
    filter <- "none"
  } else {
    dom <- "Blfrtip"
    paging <- TRUE
    filter <- "top"
  }
  
  integer_idx <- as.integer(which(sapply(df, class) == "integer"))
  
  dt <- datatable(
    df,
    filter = filter,
    class = "compact row-border stripe hover nowrap",
    extensions = "Buttons",
    caption = htmltools::tags$caption(
      style = 'caption-side: bottom; text-align: center;', caption
    ),
    options = list(
      scrollX = TRUE,
      paging = paging,
      pageLength = pageLength,
      autoWidth = TRUE,
      dom = dom,
      buttons = c("copy", "csv", "excel"),
      columnDefs = list(list(className = 'dt-center', targets = integer_idx))
    )
  )
  
  ## Numeric columns
  if (numr_cols[1] == "auto") {
    numr_cols <- names(df)[which(sapply(df, class) == "numeric")]
  }
  if (!is.null(numr_cols) & length(numr_cols) > 0) {
    dt <- dt %>% formatSignif(numr_cols, digits = 3)
  }
  
  return(dt)
}

## Misc functions
f_sci <- function(x) format(x, scientific = TRUE, digits = 2)
f_dec <- function(x) format(x, scientific = FALSE, digits = 2)


## Heatmap plot showing abundances
plot_heatmap <- function(IDs,
                         count_mat = NULL, meta_df = NULL,
                         groups = c("tissue"),
                         show_rownames = TRUE, show_colnames = FALSE,
                         cluster_rows = TRUE,
                         id_labsize = 10, ...) {
  
  ## Select groups and IDs
  fmeta <- meta_df[, groups, drop = FALSE]
  fcount_mat <- count_mat[match(IDs, rownames(count_mat)),
                          match(rownames(fmeta), colnames(count_mat))]
  fcount_mat <- as.matrix(fcount_mat)
  
  ## Arrange metadata according to the columns with included factors
  if (length(groups) == 1) fmeta <- fmeta %>% arrange(.data[[groups[1]]])
  if (length(groups) == 2) fmeta <- fmeta %>% arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]])
  if (length(groups) == 3) fmeta <- fmeta %>% arrange(.data[[groups[1]]],
                                                      .data[[groups[2]]],
                                                      .data[[groups[3]]])
  
  ## If few features are included, reduce the cell (row) height
  cellheight <- ifelse(length(IDs) > 20, NA, 20)
  id_labsize <- ifelse(length(IDs) > 40, 6, id_labsize)
  
  ## Truncate long taxon names
  row.names(fcount_mat) <- str_trunc(row.names(fcount_mat),
                                     width = 20, ellipsis = "")
  
  ## Function to create the plot
  pheatmap(fcount_mat, annotation_col = fmeta,
           cluster_rows = cluster_rows, cluster_cols = FALSE,
           show_rownames = show_rownames, show_colnames = show_colnames,
           cellheight = cellheight,
           fontsize = 9, fontsize_row = id_labsize, cex = 1,
           ...)
}