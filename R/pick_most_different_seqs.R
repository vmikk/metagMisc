
#' Pick the N most different sequences from a FASTA / XStringSet / alignment
#'
#' @param x Path to FASTA file, a DNAStringSet, or a DNAMultipleAlignment.
#' @param top_n How many sequences to select (default 3).
#' @param method "maximin" (default) or "maxsum" greedy selection.
#' @param correction Substitution model for DECIPHER::DistanceMatrix (default "Jukes-Cantor").
#' @param processors Number of CPU cores for DECIPHER (NULL = auto).
#' @param verbose Logical; pass through to DECIPHER functions.
#' @return A list with:
#'   - names: character vector of selected sequence names (in chosen order)
#'   - indices: integer indices into the original input order (if available)
#'   - subset: DNAStringSet of selected sequences (original, unaligned)
#'   - dist: numeric matrix of pairwise distances used for selection
#' @examples
#' # res <- pick_most_different("seqs.fasta", top_n = 5)
#' # res$subset  # DNAStringSet of the 5 most different sequences

## Pick the N most different sequences
pick_most_different_seqs <- function(x, top_n = 3,
  method = c("maximin", "maxsum"),
  correction = "JC69", processors = NULL, verbose = FALSE) {

  method <- match.arg(method)

  # x = DNAStringSet
  # method = "maximin" (default, farthest-point sampling) or "maxsum" (greedy selection, "most distant-on-average")

  ss <- x
  if (!inherits(ss, "DNAStringSet")) stop("x must be a DNAStringSet.")
  src_names <- names(ss)

  if (length(ss) < 2L) stop("Need at least two sequences.")
  if (is.null(src_names)) src_names <- paste0("seq", seq_along(ss))

  ## Multiple sequence alignment
  if (verbose) message("Aligning ", length(ss), " sequences ...")
  aln <- DECIPHER::AlignSeqs(
      myXStringSet = ss,
      processors = processors,
      verbose = verbose)   # DNAStringSet with gaps

  ## Distance matrix
  if (verbose) message("Computing distance matrix ...")
  dd <- DECIPHER::DistanceMatrix(
    myXStringSet        = aln,
    method              = "overlap",
    type                = "matrix",
    includeTerminalGaps      = FALSE,
    penalizeGapLetterMatches = FALSE,  # ignore gaps paired with letters
    correction          = correction,
    processors          = processors,
    verbose             = verbose)

  dd <- as.matrix(dd)
  diag(dd) <- 0

  ## Handle NA distances (e.g., no overlap)
  ## TO DO - drop sequences with too many NAs ???
  if (all(is.na(dd))) {
    stop("All pairwise distances are NA (no overlap among sequences). Try different alignment settings.")
  }
  if (anyNA(dd)) {
    if (verbose) message("NA distances found; treating them as the current maximum finite distance.")
    max_finite <- max(dd, na.rm = TRUE)
    dd[ is.na(dd) ] <- max_finite
  }

  n <- nrow(dd)
  rn <- rownames(dd)
  if (is.null(rn)) rn <- src_names

  ## If asking for everything, return as-is
  if (n <= top_n) {
    res <- list(
      names = rn,
      indices = seq_len(n),
      subset = ss,
      dist = dd)

    return(res)
  }

  ## Special case: top_n == 1 -> pick the sequence most distant on average
  if (top_n == 1L) {
    i <- which.max(rowMeans(dd))
    res <- list(
      names   = rn[i],
      indices = i,
      subset  = ss[i],
      dist    = matrix(0, 1, 1, dimnames = list(rn[i], rn[i])) )
    return(res)
  }

  ## Seed with the farthest pair from the upper triangle
  upper <- dd
  upper[lower.tri(upper, diag = TRUE)] <- -Inf
  ij <- which(upper == max(upper, na.rm = TRUE), arr.ind = TRUE)[1L, ]
  selected <- c(ij[1L], ij[2L])
  names(selected) <- rn[ selected ]

  ## Incremental trackers
  cur_min <- pmin(dd[, selected[1L]], dd[, selected[2L]])   # for maximin
  cur_sum <- dd[, selected[1L]] + dd[, selected[2L]]        # for maxsum
  cur_min[selected] <- -Inf
  cur_sum[selected] <- -Inf

  ## Greedy expansion
  while (length(selected) < top_n) {
    pick <- if (method == "maximin") which.max(cur_min) else which.max(cur_sum)
    selected <- c(selected, pick)

    ## Update trackers
    colp <- dd[, pick]
    cur_min <- pmin(cur_min, colp)
    cur_sum <- cur_sum + colp

    ## Exclude chosen from future picks
    cur_min[selected] <- -Inf
    cur_sum[selected] <- -Inf
  }

  res <- list(
    names   = rn[ selected ],
    indices = selected,
    subset  = ss[selected],                       # original (unaligned) sequences
    dist    = dd[selected, selected, drop = FALSE] )

  return(res)
}


