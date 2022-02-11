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
  if (numr_cols == "auto") {
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