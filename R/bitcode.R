
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

