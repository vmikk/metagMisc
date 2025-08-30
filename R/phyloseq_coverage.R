
#' @title Estimate the observed abundance-based sample coverage for phyloseq object
#' @description phyloseq_coverage estimates the sample completeness for the individual-based
#' abundance data (number of sequencing reads) stored in 'phyloseq'-class objects.
#' @param physeq A phyloseq-class object
#' @param correct_singletons Logical; if TRUE, singleton counts will be corrected with modified Good–Turing frequency formula (Chiu, Chao 2016)
#' @param add_attr Logical; if TRUE, additional attributes (list of species abundances and singleton correction flag) will be added to the results
#' @details Coverage represents a measure of sample completeness and is defined as the proportion of
#' the total number of individuals in a community that belong to the species represented in the sample.
#' Coverage complement (1 - Coverage) gives the proportion of the community belonging to unsampled
#' species or the "coverage deficit" (Chao, Jost, 2012).
#' Estimation of coverage is based on the number of singletons and doubletons in the sample.
#' @return Data frame with coverage estimates for each sample
#' @export
#' @references
#' Chao A, Jost L. (2012) Coverage-based rarefaction and extrapolation: standardizing samples by completeness rather than size // Ecology 93(12): 2533–2547. DOI: 10.1890/11-1952.1
#' @examples
#' data("esophagus")
#' phyloseq_coverage(esophagus)
#'
phyloseq_coverage <- function(physeq, correct_singletons = FALSE, add_attr = T){

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        phyloseq_otu_to_df(physeq, taxa_as_rows = T),
        correct_singletons = correct_singletons)

  ## Estimate sample coverages
  res <- plyr::ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
  colnames(res) <- c("SampleID", "SampleCoverage")

  ## Add attributes
  if(add_attr){
    attr(res, "x") <- x
    attr(res, "correct_singletons") <- correct_singletons
  }

  return(res)
}



## Estimate the required sample size for a particular coverage
## the code is based on iNEXT:::invChat.Ind by Johnson Hsieh (d76e3b8, Nov 12, 2016)
# https://github.com/JohnsonHsieh/iNEXT/blob/de46aeacb4433c539b2880df376e87b44bc1723c/R/invChat.R#L1
coverage_to_samplesize <- function(x, coverage = 0.95, add_attr = F){
  # x = vector of species abundances
  # coverage = the desired sample completness that we want to achieve

  # iNEXT:::invChat.Ind(x, C = coverage)$m

  ## Total number of reads and the observed sample coverage
  n <- sum(x)
  refC <- iNEXT:::Chat.Ind(x, n)

  ## Interpolation
  f <- function(m, C) abs(iNEXT:::Chat.Ind(x, m) - C)
  if(refC > coverage){
    opt <- optimize(f, C = coverage, lower = 0, upper = sum(x))
    mm <- opt$minimum
    mm <- round(mm)
  }

  ## Extrapolation
  if(refC <= coverage){
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if(f1>0 & f2>0)  {A <- (n-1)*f1/((n-1)*f1+2*f2)}
    if(f1>1 & f2==0) {A <- (n-1)*(f1-1)/((n-1)*(f1-1)+2)}
    if(f1==1 & f2==0){A <- 1}
    if(f1==0 & f2==0){A <- 1}
    mm <- (log(n/f1)+log(1-coverage))/log(A)-1
    mm <- n + mm
    mm <- round(mm)
  }

  ## Add attributes
  if(add_attr){
    if(refC > coverage) { attr(mm, "method") <- "interpolated" }
    if(refC <= coverage){ attr(mm, "method") <- "extrapolated" }
    attr(mm, "ObservedCoverage") <- refC
    attr(mm, "RequestedCoverage") <- coverage
  }

  return(mm)
}
## Example:
# abunds <- c(48,21,16,15,14,6,6,2,2,2,1,1,1,1,1)
# coverage_to_samplesize(abunds, coverage = 0.9, add_attr = T)   # interpolated sample size
# coverage_to_samplesize(abunds, coverage = 0.97, add_attr = T)  # extraplation




#' @title Coverage-based rarefaction
#' @description This function performs coverage-based rarefaction (interpolation) based on the analytical approach proposed by Chao and Jost (2012).
#' @param physeq A phyloseq-class object
#' @param coverage Numeric value for a particular sample coverage (between 0 and 1)
#' @param iter Number of rarefication iterations
#' @param replace Logical, whether to sample with replacement (TRUE) or without replacement (FALSE, default)
#' @param correct_singletons Logical; if TRUE, singleton counts will be corrected with modified Good–Turing frequency formula (Chiu, Chao 2016)
#' @param seeds Integer vector used for the reproducible random subsampling (should be of the same length as the number of iterations)
#' @param multithread Logical or integer; if TRUE, attempts to run the function on multiple cores; integer defines the number of cores to use (if it is set to TRUE, all cores will be used)
#' @param drop_lowcoverage Logical; if TRUE, samples with coverage lower than selected value will be removed (default, FALSE)
#' @param ... Additional arguments will be passed to \code{\link[phyloseq]{rarefy_even_depth}}
#' @details
#' Samples standardized by size will have different degrees of completness.
#' When we compare samples with the same coverage, we are making sure that samples are equally complete
#' and that the unsampled species constitute the same proportion of the total individuals in each community (Chao, Jost, 2012).
#'
#' @return List of rarefied phyloseq-objects (or a single phyloseq object if iter = 1)
#' @export
#'
#' @examples
#' # Load data
#' data("esophagus")
#' 
#' # Coverage-based rarefaction
#' eso_raref <- phyloseq_coverage_raref(physeq = esophagus, iter = 1, coverage = 0.8)
#' 
#' # Perform coverage-based rarefaction multiple times (iter = 5)
#' eso_raref2 <- phyloseq_coverage_raref(physeq = esophagus, iter = 5, coverage = 0.8)
#' 
phyloseq_coverage_raref <- function(physeq, coverage = NULL, iter = 1, replace = F,
  correct_singletons = FALSE, seeds = NULL, multithread = F, drop_lowcoverage = F, ...){

  ## Prepare seed values for rarefaction
  if(is.null(seeds)){ seeds <- 1:iter }

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        phyloseq_otu_to_df(physeq, taxa_as_rows = T),
        correct_singletons = correct_singletons)

  ## Estimate the observed sample coverages
  SC <- plyr::ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
  colnames(SC) <- c("SampleID", "SampleCoverage")
  SC$SampleID <- as.character(SC$SampleID)

  ## Select the lowest observed coverage
  if(is.null(coverage)){
    coverage <- min(SC$SampleCoverage)
    cat("Coverage value was set to the minimum observed value across all samples (", coverage, ")\n", sep = "")
  }

  ## Data validation
  if(any(SC$SampleCoverage < coverage)){
    if(drop_lowcoverage == FALSE){
      stop("There are not enough data to reach the required coverage for some samples.\n")
    }

    if(drop_lowcoverage == TRUE){
      lowcov <- SC$SampleCoverage < coverage
      nlow <- sum(lowcov)
      warning("Samples with coverage lower than the selected threshold were discared (n = ", nlow, ").\n")

      ## Remove samples
      x <- x[ !lowcov ]
      physeq <- phyloseq::prune_samples(SC$SampleID[ !lowcov ], physeq)
    }
  } # End of data validation

  ## Estimate the required sample sizes
  RSZ <- plyr::ldply(.data = x, .fun = coverage_to_samplesize, coverage = coverage, add_attr = F)
  colnames(RSZ) <- c("SampleID", "RequiredSize")
  rownames(RSZ) <- RSZ$SampleID
  RSZ$SampleID <- as.character(RSZ$SampleID)

  ## Split the samples
  phys <- phyloseq_sep_samp(physeq, drop_zeroes = FALSE)

  ##### Rarefy each sample

  ## Single rarefaction
  if(iter == 1){

    res <- plyr::mlply(
      .data = RSZ$SampleID,       # input = sample names
      .fun = function(z, ...){

        ## Extract phyloseq object for a single sample
        pp <- phys[[z]]

        ## Extract the required sample size
        SampSize <- RSZ[z, "RequiredSize"]

        ## Rarefaction
        rz <- phyloseq::rarefy_even_depth(pp, sample.size=SampSize, verbose = F, trimOTUs = FALSE, rngseed = seeds, ...)
        return(rz)
      },
      replace=replace,
      ...)  # pass additional arguments to  phyloseq::rarefy_even_depth


      ## Combine rarefied samples into a single phyloseq object
      RES <- res[[1]]
      for(i in 2:length(res)){
        RES <- phyloseq::merge_phyloseq(RES, res[[i]])
      }
      rm(i)

  } # end of iter == 1


  ## Multiple rarefaction
  if(iter > 1){

    ############### Cluster setup

    ## Progress bar type for single-threaded plyr functions
    progr <- "text"
    parall <- FALSE

    ## Check if foreach and doParallel packages are available
    if(multithread){
      if(!requireNamespace("foreach", quietly = TRUE)){ stop("foreach package required for parallel plyr operation.\n") }
      if(!requireNamespace("doParallel", quietly = TRUE)){ stop("doParallel package is required.\n") }
    }

    ## Check the platform type and specify the number of cores to use
    if(multithread && .Platform$OS.type == "unix"){
      ncores <- parallel::detectCores()
      if(is.numeric(multithread)){ ncores <- multithread }
      if(is.na(ncores) | is.null(ncores)){ ncores <- 1 }
    } else {
      ncores <- 1
      if(multithread && .Platform$OS.type=="windows"){
        warning("Multithreading has been DISABLED, as forking is not supported on .Platform$OS.type 'windows'.\n")
      }
    }

    ## Setup cluster
    if(ncores > 1){

      ## plyr arguments for parallel execution
      progr <- "none"
      parall <- TRUE

      ## Disable load balancing
      paropts <- list(preschedule=TRUE)

      ## Start the cluster
      cl <- parallel::makeCluster(ncores)

      ## Register the parallel backend
      doParallel::registerDoParallel(cl)

      ## Load packages on cluster nodes
      parallel::clusterEvalQ(cl, library("phyloseq"))

      ## Send useful objects to the workers
      vars <- c("RSZ", "replace", "phys")
      parallel::clusterExport(cl=cl, varlist=vars, envir = environment())
    }

    ############### End of cluster setup

    ## SampleID - Seed pairs
    SSP <- expand.grid(rngseed = seeds, SampleID = RSZ$SampleID, stringsAsFactors = F)

    ## Rarefy
    res <- plyr::mlply(
      .data = SSP,
      .fun = function(rngseed = rngseed, SampleID = SampleID, ...){

          ## Extract phyloseq object for a single sample
          pp <- phys[[ SampleID ]]

          ## Extract the required sample size
          SampSize <- RSZ[SampleID, "RequiredSize"]

          ## Rarefaction
          rz <- phyloseq::rarefy_even_depth(pp, sample.size=SampSize, verbose = F, trimOTUs = FALSE, rngseed = rngseed, ...)
          return(rz)
        },
      .progress = progr,
      .parallel = parall,
      replace=replace,
      ...)  # pass additional arguments to  phyloseq::rarefy_even_depth


    ## Combine samples form the same iteration into a single list
    res <- split(x = res, f = SSP$rngseed)

    ## Merge samples (within the same iteration)
    RES <- plyr::llply(.data = res, .fun = function(z){
      tmp <- z[[1]]
      for(i in 2:length(z)){
        tmp <- phyloseq::merge_phyloseq(tmp, z[[i]])
      }
      return(tmp)
    })


    ## Stop the cluster
    if(ncores > 1){
      parallel::stopCluster(cl)
      rm(cl)
    }

  } # end of iter > 1


  ## Add attributes
  attr(RES, "SampleCoverage") <- RSZ

  return(RES)
}

