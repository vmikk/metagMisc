
# Estimate the observed abundance-based sample coverage for phyloseq object
phyloseq_coverage <- function(physeq, correct_singletons = FALSE, add_attr = T){
  
  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate sample coverages
  res <- ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
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
# Example:  coverage_to_samplesize(x, coverage = 0.9, add_attr = T)



## Perform rarefaction with 
phyloseq_coverage_raref <- function(physeq, coverage = 0.95, iter = 1, replace = F, 
  correct_singletons = FALSE, seeds = NULL, multithread = F, ...){
  # ... passed to rarefy_even_depth

  ## Prepare seed values for rarefaction
  if(is.null(seeds)){ seeds <- 1:iter }

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate the observed sample coverages
  SC <- ldply(.data = x, .fun = function(z){ iNEXT:::Chat.Ind(z, sum(z)) })
  colnames(SC) <- c("SampleID", "SampleCoverage")

  ## Data validation
  if(any(SC$SampleCoverage < coverage)){
    stop("There are not enough data to reach the required coverage for some samples.\n")
  }

  ## Estimate the required sample sizes
  RSZ <- ldply(.data = x, .fun = coverage_to_samplesize, coverage = coverage, add_attr = F)
  colnames(RSZ) <- c("SampleID", "RequiredSize")
  rownames(RSZ) <- RSZ$SampleID
  RSZ$SampleID <- as.character(RSZ$SampleID)

  ## Split the samples
  phys <- phyloseq_sep_samp(physeq, drop_zeroes = FALSE)

  ##### Rarefy each sample
  
  ## Single rarefaction
  if(iter == 1){
  
    res <- mlply(
      .data = RSZ$SampleID,       # input = sample names
      .fun = function(z, ...){ 
  
        ## Extract phyloseq object for a single sample
        pp <- phys[[z]]
  
        ## Extract the required sample size
        SampSize <- RSZ[z, "RequiredSize"]
  
        ## Rarefaction
        rz <- rarefy_even_depth(pp, sample.size=SampSize, verbose = F, trimOTUs = FALSE, rngseed = seeds, ...)
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
    res <- mlply(
      .data = SSP,
      .fun = function(rngseed = rngseed, SampleID = SampleID, ...){
 
          ## Extract phyloseq object for a single sample
          pp <- phys[[ SampleID ]]
    
          ## Extract the required sample size
          SampSize <- RSZ[SampleID, "RequiredSize"]
    
          ## Rarefaction
          rz <- rarefy_even_depth(pp, sample.size=SampSize, verbose = F, trimOTUs = FALSE, rngseed = rngseed, ...)
          return(rz)
        },
      .progress = progr,
      .parallel = parall,
      replace=replace,
      ...)  # pass additional arguments to  phyloseq::rarefy_even_depth


    ## Combine samples form the same iteration into a single list
    res <- split(x = res, f = SSP$rngseed)

    ## Merge samples (within the same iteration)
    RES <- llply(.data = res, .fun = function(z){
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

