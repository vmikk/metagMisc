
#' @title RTK-based rarefaction
#'
#' @param physeq A phyloseq-class object
#' @param SampSize Rarefaction depth (number of reads to sample)
#' @param MinSizeTreshold Remove samples with number of reads less then this treshold
#' @param iter Number of rarefication iterations
#' @param trimOTUs Logical, if TRUE (default), OTUs that have a count of zero in every sample will be removed
#' @param ... Additional arguments will be passed to \code{\link[rtk]{rtk}}
#'
#' @return List of rarefied phyloseq-objects.
#' @export
#'
#' @examples
#'
raref_rtk <- function(physeq, SampSize = NULL, MinSizeTreshold = NULL,
  iter = 1000, trimOTUs = TRUE, ...){

  ## Filter samples by number of reads
  if(!is.null(MinSizeTreshold)){ x <- phyloseq::prune_samples(phyloseq::sample_sums(x) >= MinSizeTreshold, x) }

  ## Define rarefication depth
  if(is.null(SampSize)){ SampSize <- round( 0.9*min(phyloseq::sample_sums(x)) ) }

  ## Extract OTU abundance table (should be a matrix)
  abunds <- as( phyloseq::otu_table(physeq), "matrix" )

  ## Define margins
  TR <- phyloseq::taxa_are_rows(physeq)
  if(TR == TRUE) { MAR <- 2 }
  if(TR == FALSE){ MAR <- 1 }
  #   columns represent single samples (margin=2)
  #   rows are assumed to be samples (margin=1)

  ## Perform rarefaction
  rar <- rtk_mod(abunds, repeats = iter, depth = SampSize,
    margin = MAR, verbose = FALSE, threads = 1, ...)

  ## Function to replace OTU table in phyloseq object
  subst_otu_tab <- function(phys, newotu, taxrows, drop_zeros = FALSE){
    phyloseq::otu_table(phys) <- phyloseq::otu_table(newotu, taxa_are_rows = taxrows)

    ## Remove OTUs from the dataset that are no longer observed in any sample
    if(drop_zeros == TRUE){
      phys <- phyloseq::prune_taxa(phyloseq::taxa_sums(phys) > 0, phys)
    }

    return(phys)
  }
  ## e.g. subst_otu_tab(physeq, newotu = rar[[1]])

  ## Re-create phyloseq object (single depth case)
  if(length(SampSize) == 1){
    res <- llply(
      .data = rar,
      .fun = function(z){
        subst_otu_tab(phys = physeq, newotu = z, taxrows = TR, drop_zeros = trimOTUs)
      })
  }

  ## Re-create phyloseq object (multiple depth case)
  # if(length(SampSize) > 1){
  #   ///////////////////       - TO DO
  # }


  ## Add rarefaction parameters as attributes to the phyloseq object
  attr(res, which = "RarefactionDepth") <- SampSize
  attr(res, which = "taxa_are_rows") <- TR

  return(res)
}


## Modified RTK function that returns only rarefied matrices (without diversity estimation)
# based on https://github.com/hildebra/Rarefaction/blob/05b8f37e27d7c023334846df8e690dbbd833a4a9/r-package/rtk/R/rarefaction.R#L13
# Original author - Paul Saary
rtk_mod <- function(input, repeats = 10, depth = 0, margin = 2,
    verbose = FALSE, threads = 1, lowmem = FALSE){   # tmpdir = NULL

  ## Pass 1:x to Cpp as colnames
  removeCnames <- FALSE
  removeRnames <- FALSE

  ## Return all rarefied matrices
  ReturnMatrix <- repeats

  ## Sort depths
  depth <- sort(as.numeric(depth))

  ## Convert dataframe to matrix
  # if(class(input) == "data.frame"){
  #   input <- as.matrix(input)
  # }

  ## Validate that the matrix is numeric
  if(!is.numeric(input)){
    stop("The supplied matrix object is not numeric. Please check your input matrix.")
  }
  if(any(is.na(input))){
    stop("The input data contains NA values. Please sanitize your input first.")
  }

  if(is.null(colnames(input))){
    colnames(input) <- paste("col ", seq(1:ncol(input)), sep="")
    removeCnames <- TRUE
  }
  if(is.null(rownames(input))){
    rownames(input) <- paste("row ", seq(1:nrow(input)), sep="")
    removeRnames <- TRUE
  }

  ## Call the actual software
  result <- rtk:::rcpp_rarefaction(
    input = "",
    rMatrix = input,
    inColNames = colnames(input),
    inRowNames = rownames(input),
    repeats = repeats,
    depth = depth,
    NoOfMatrices = ReturnMatrix,
    verbose = verbose,
    threads = threads,
    margin = margin,
    tmpDir = "NULL",
    lowmem = lowmem)

  ## Call the garbage collecor
  gc()


  result <- lapply(result, function(res){

    ## Remove names, if there werent any
    if(removeRnames == TRUE && removeCnames == TRUE){
      res$raremat <- lapply(res$raremat, unname)
    } else {
      if(removeRnames == TRUE){
        res$raremat <- lapply(res$raremat, function(x){ rownames(x) <- NULL; return(x) })
      }
      if(removeCnames == TRUE){
        res$raremat <- lapply(res$raremat, function(x){ colnames(x) <- NULL; return(x) })
      }
    }

    if(length(res$skipped) > 0){
      warning(paste(length(res$skipped), "samples where skipped because the depth was greater than the number of elements in the sample."))
    }

    ## Remove everything except rarefied matrices
    res <- res$raremat

    return(res)
  })


  if(length(depth) == 1){
    result <- result[[1]]        # only 1 sampling depth
  } else {
    names(result) <- depth       # multiple sampling depths
  }

  attr(result, which = "depths") <- depth
  attr(result, which = "repeats") <- repeats

  gc()
  return(result)
}
