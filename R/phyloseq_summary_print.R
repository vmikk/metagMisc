
## Helper functions to print phyloseq summary in a human-readable format

#' @noRd
phyloseq_summary_parameter_type <- function(parameter) {
  count_parameters <- c(
    "Number of samples",
    "Number of OTUs",
    "Total number of reads",
    "Number of singletons",
    "Data sparsity (number of zeros)",
    "Min total OTU abundance",
    "Max total OTU abundance",
    "Min total sample abundance",
    "Max total sample abundance"
  )

  is_percent <- grepl("percentage", parameter, ignore.case = TRUE) |
    grepl("occurrence, percents$", parameter)

  out <- rep("continuous", length(parameter))
  out[parameter %in% count_parameters] <- "count"
  out[is_percent] <- "percent"
  out
}

#' @noRd
phyloseq_summary_format_number <- function(x, digits = 4L) {
  out <- rep(NA_character_, length(x))

  finite <- is.finite(x)
  if (!any(finite)) {
    out[is.na(x)] <- NA_character_
    out[is.nan(x)] <- "NaN"
    out[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, "Inf", "-Inf")
    return(out)
  }

  rounded <- round(x[finite])
  integerish <- abs(x[finite] - rounded) < sqrt(.Machine$double.eps)
  finite_idx <- which(finite)

  if (any(integerish)) {
    out[finite_idx[integerish]] <- prettyNum(
      rounded[integerish],
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  }

  if (any(!integerish)) {
    vals <- signif(x[finite][!integerish], digits = digits)
    out[finite_idx[!integerish]] <- prettyNum(
      vals,
      digits = digits,
      big.mark = ",",
      scientific = FALSE,
      preserve.width = "none",
      drop0trailing = TRUE,
      trim = TRUE
    )
  }

  out[is.nan(x)] <- "NaN"
  out[is.infinite(x)] <- ifelse(x[is.infinite(x)] > 0, "Inf", "-Inf")
  out
}

#' @noRd
phyloseq_summary_format_column <- function(values, parameter_type) {
  out <- rep(NA_character_, length(values))

  is_count <- parameter_type == "count"
  if (any(is_count)) {
    out[is_count] <- prettyNum(
      round(values[is_count]),
      big.mark = ",",
      scientific = FALSE,
      trim = TRUE
    )
  }

  is_other <- !is_count
  if (any(is_other)) {
    out[is_other] <- phyloseq_summary_format_number(values[is_other], digits = 4L)
  }

  is_percent <- parameter_type == "percent"
  out[is_percent] <- paste0(out[is_percent], "%")
  out
}

#' @noRd
phyloseq_summary_wrap_lines <- function(text, width) {
  if (is.na(text) || is.null(text)) {
    return(NA_character_)
  }

  wrapped <- strwrap(text, width = width, simplify = FALSE)
  if (length(wrapped) == 0L || length(wrapped[[1]]) == 0L) {
    return("")
  }

  wrapped[[1]]
}

#' @noRd
phyloseq_summary_style_fns <- function(use_color = TRUE) {
  identity_fn <- function(x) x

  if (!isTRUE(use_color) ||
      !requireNamespace("crayon", quietly = TRUE) ||
      !isTRUE(crayon::has_color())) {
    return(list(
      header     = identity_fn,
      parameter  = identity_fn,
      count      = identity_fn,
      percent    = identity_fn,
      continuous = identity_fn
    ))
  }

  list(
    header     = crayon::blue$bold,
    parameter  = crayon::silver,
    count      = crayon::green,
    percent    = crayon::yellow,
    continuous = identity_fn
  )
}

#' @noRd
format.phyloseq_summary_wide <- function(x, ..., use_color = TRUE) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  stopifnot("Parameter" %in% colnames(df))

  parameter_type <- phyloseq_summary_parameter_type(df$Parameter)
  value_cols <- setdiff(colnames(df), "Parameter")

  display <- data.frame(Parameter = df$Parameter, stringsAsFactors = FALSE)
  for (col in value_cols) {
    display[[col]] <- phyloseq_summary_format_column(df[[col]], parameter_type)
  }

  plain_display <- display
  sep <- "  "

  value_widths <- vapply(colnames(plain_display)[-1], function(col) {
    max(nchar(c(col, plain_display[[col]]), type = "width"), na.rm = TRUE)
  }, integer(1))

  total_value_width <- 0L
  if (length(value_widths) > 0L) {
    total_value_width <- sum(value_widths) + (length(value_widths) * nchar(sep, type = "width"))
  }

  table_width <- getOption("width", 80L)
  parameter_width <- max(
    nchar("Parameter", type = "width"),
    min(
      max(nchar(plain_display$Parameter, type = "width"), na.rm = TRUE),
      max(20L, table_width - total_value_width)
    )
  )

  styles <- phyloseq_summary_style_fns(use_color = use_color)

  header_cells <- c(
    format("Parameter", width = parameter_width, justify = "left"),
    vapply(seq_along(value_cols), function(idx) {
      format(value_cols[[idx]], width = value_widths[[idx]], justify = "right")
    }, character(1))
  )
  header_cells[1] <- styles$header(header_cells[1])
  if (length(value_cols) > 0L) {
    header_cells[-1] <- vapply(header_cells[-1], styles$header, character(1))
  }

  lines <- paste(header_cells, collapse = sep)

  for (row_idx in seq_len(nrow(plain_display))) {
    wrapped_parameter <- phyloseq_summary_wrap_lines(
      plain_display$Parameter[[row_idx]],
      width = parameter_width
    )

    numeric_cells <- vapply(seq_along(value_cols), function(idx) {
      format(
        plain_display[row_idx, value_cols[[idx]], drop = TRUE],
        width = value_widths[[idx]],
        justify = "right"
      )
    }, character(1))

    numeric_style <- styles[[parameter_type[[row_idx]]]]
    if (length(numeric_cells) > 0L) {
      numeric_cells <- vapply(numeric_cells, numeric_style, character(1))
    }

    for (line_idx in seq_along(wrapped_parameter)) {
      parameter_cell <- format(
        wrapped_parameter[[line_idx]],
        width = parameter_width,
        justify = "left"
      )
      parameter_cell <- styles$parameter(parameter_cell)

      if (line_idx == 1L) {
        row_cells <- c(parameter_cell, numeric_cells)
      } else {
        blank_cells <- vapply(value_widths, function(width) {
          format("", width = width, justify = "right")
        }, character(1))
        row_cells <- c(parameter_cell, blank_cells)
      }

      lines <- c(lines, paste(row_cells, collapse = sep))
    }
  }

  lines
}

#' @noRd
print.phyloseq_summary_wide <- function(x, ..., use_color = TRUE) {
  cat(format(x, ..., use_color = use_color), sep = "\n")
  cat("\n")
  invisible(x)
}
