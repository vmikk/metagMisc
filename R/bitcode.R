
#' Encode combinations of discrete variables into integer codes
#'
#' @description
#' Encodes joint states of 2+ discrete (categorical) variables into a single
#' integer using a mixed-radix representation. This is useful for building
#' compact IDs for combinations (e.g., taxonomy rank tuples, multi-factor groups),
#' and is closely related to enumerating states of a product alphabet in
#' information theory.
#'
#' @param ... Vectors (all the same length) representing discrete variables to
#'   encode. Alternatively, provide \code{x}.
#' @param x A \code{data.frame}, \code{matrix}, or \code{list} of equal-length
#'   vectors to encode. If provided, \code{...} must be empty.
#' @param levels Optional \code{list} defining the allowed levels (and their
#'   order) for each variable. If unnamed, must have the same length and order
#'   as the input variables. If named, names must match the variable names.
#'   By default, factors use \code{levels(x)}, while other vectors use
#'   \code{unique(x)} in order of appearance (excluding \code{NA}).
#' @param start Integer scalar; the first code value. Use \code{0} (default) to
#'   get codes in \code{0, 1, ..., prod(bases) - 1}; use \code{1} for 1-based IDs.
#' @param order Character scalar; either \code{"most_significant_first"} (default)
#'   or \code{"least_significant_first"}. With \code{"most_significant_first"},
#'   the first variable is treated as the most significant "digit":
#'   \eqn{code = i_1 \prod_{j=2}^k b_j + i_2 \prod_{j=3}^k b_j + \cdots + i_k}.
#' @param na_code Value to use for rows containing missing values in any input
#'   variable. Default is \code{NA} (a missing code).
#' @param output Character scalar; \code{"auto"} (default) returns integer when
#'   safe, otherwise numeric. Use \code{"integer"} or \code{"numeric"} to force.
#'
#' @return A vector of codes (integer or numeric) of the same length as the input
#'   variables, with attributes describing the encoding:
#'   \itemize{
#'     \item \code{bitcode_levels}: list of level sets per variable
#'     \item \code{bitcode_bases}: integer vector of base sizes per variable
#'     \item \code{bitcode_order}: encoding order
#'     \item \code{bitcode_start}: code start offset
#'     \item \code{bitcode_is_factor}: logical vector (whether each input was a factor)
#'     \item \code{bitcode_varnames}: variable names used for decoding
#'   }
#'
#' @details
#' This is a mixed-radix (positional) encoding of a Cartesian product of
#' alphabets, i.e. a reversible map between tuples and integers as long as the
#' level sets are fixed. By default, level sets are derived from the data, which
#' makes codes stable within that dataset; for stable codes across datasets, pass
#' \code{levels}.
#'
#' @examples
#' # Two variables (most significant first)
#' a <- c("A", "A", "B", "B")
#' b <- c("x", "y", "x", "y")
#' code <- bitcode(a, b, start = 0)
#' code
#'
#' # Decode back (uses attributes stored on 'code')
#' bitdecode(code)
#'
#' # Factors: uses factor levels by default
#' f <- factor(c("low", "high", "low"), levels = c("low", "med", "high"))
#' g <- factor(c("x", "y", "x"))
#' code2 <- bitcode(f, g, start = 1)
#' bitdecode(code2)
#'
#' @export
bitcode <- function(
    ...,       # vectors to encode
    x = NULL,
    levels = NULL,
    start = 0L,
    order = c("most_significant_first", "least_significant_first"),
    na_code = NA,
    output = c("auto", "integer", "numeric")) {

  order <- match.arg(order)
  output <- match.arg(output)

  if (!is.null(x) && length(list(...)) > 0) {
    stop("Provide either 'x' or '...', not both.")
  }

  vars <- list(...)
  if (!is.null(x)) {
    if (is.data.frame(x)) {
      vars <- as.list(x)
      if (is.null(names(vars))) names(vars) <- colnames(x)
    } else if (is.matrix(x)) {
      vars <- as.data.frame(x, stringsAsFactors = FALSE)
      vars <- as.list(vars)
      if (is.null(names(vars))) names(vars) <- colnames(x)
    } else if (is.list(x)) {
      vars <- x
    } else {
      stop("'x' must be a data.frame, matrix, or list.")
    }
  }

  if (length(vars) < 2) stop("Provide at least two variables to encode.")

  n <- length(vars[[1]])
  if (!all(vapply(vars, length, integer(1)) == n)) {
    stop("All input variables must have the same length.")
  }

  varnames <- names(vars)
  if (is.null(varnames) || any(varnames == "")) {
    varnames <- paste0("V", seq_along(vars))
  }

  is_factor <- vapply(vars, is.factor, logical(1))

  ## Build / validate level sets
  lvl <- vector("list", length(vars))
  names(lvl) <- varnames

  if (!is.null(levels)) {
    if (!is.list(levels)) stop("'levels' must be a list or NULL.")
    if (!is.null(names(levels)) && any(names(levels) != "")) {
      missing_names <- setdiff(varnames, names(levels))
      if (length(missing_names) > 0) {
        stop("Named 'levels' is missing entries for: ", paste(missing_names, collapse = ", "))
      }
      for (nm in varnames) lvl[[nm]] <- levels[[nm]]
    } else {
      if (length(levels) != length(vars)) {
        stop("Unnamed 'levels' must have the same length as the number of variables.")
      }
      for (j in seq_along(vars)) lvl[[j]] <- levels[[j]]
      names(lvl) <- varnames
    }
  } else {
    for (j in seq_along(vars)) {
      v <- vars[[j]]
      if (is.factor(v)) {
        lvl[[j]] <- base::levels(v)
      } else {
        u <- unique(v)
        u <- u[!is.na(u)]
        lvl[[j]] <- u
      }
    }
  }

  bases <- vapply(lvl, length, integer(1))
  if (any(bases <= 0)) stop("Each variable must have at least one level.")

  ## Map each variable to 0..(base-1)
  idx <- matrix(NA_real_, nrow = n, ncol = length(vars))
  colnames(idx) <- varnames
  for (j in seq_along(vars)) {
    v <- vars[[j]]
    lev_j <- lvl[[j]]
    m <- match(v, lev_j)
    unknown <- !is.na(v) & is.na(m)
    if (any(unknown)) {
      stop("Values not present in the provided/derived levels for variable '",
           varnames[[j]], "': ", paste(utils::head(unique(v[unknown]), 5), collapse = ", "),
           if (length(unique(v[unknown])) > 5) ", ..." else "")
    }
    idx[, j] <- ifelse(is.na(m), NA_real_, m - 1)
  }

  ## Compute weights
  k <- length(vars)
  weights <- switch(
    order,
    most_significant_first = c(rev(cumprod(rev(bases[-1]))), 1),
    least_significant_first = cumprod(c(1, bases[-k]))
  )

  ## Determine safe output type
  max_code <- prod(as.double(bases)) - 1 + as.double(start)
  can_be_integer <- is.finite(max_code) && max_code <= .Machine$integer.max &&
    isTRUE(all.equal(start, as.integer(start)))

  if (output == "integer" && !can_be_integer) {
    warning("Requested output='integer' but codes may overflow integer; returning numeric instead.")
    output <- "numeric"
  }
  if (output == "auto") output <- if (can_be_integer) "integer" else "numeric"

  out <- rep(na_code, n)
  ok <- !apply(is.na(idx), 1, any)
  if (any(ok)) {
    z <- idx[ok, , drop = FALSE] %*% matrix(weights, ncol = 1)
    z <- as.vector(z) + as.double(start)
    out[ok] <- if (output == "integer") as.integer(z) else as.double(z)
  }

  attr(out, "bitcode_levels") <- lvl
  attr(out, "bitcode_bases") <- bases
  attr(out, "bitcode_order") <- order
  attr(out, "bitcode_start") <- as.integer(start)
  attr(out, "bitcode_is_factor") <- is_factor
  attr(out, "bitcode_varnames") <- varnames
  class(out) <- unique(c("bitcode", class(out)))
  out
}

