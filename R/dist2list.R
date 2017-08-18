
## Convert distance matrix to pairwised list
dist2list <- function (dist, tri=TRUE) {
  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }

  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)

  if(tri == TRUE){    # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
  }

  return(res)
}
