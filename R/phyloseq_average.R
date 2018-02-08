
#' @title Average relative OTU abundances.
#' @description This function implements OTU abundance averaging following CoDa (Compositional Data Analysis) workflow.
#' @param physeq A phyloseq-class object
#' @param avg_type Averaging type ("acomp" for Aitchison CoDa approach; "arithmetic" for simple arithmetic mean)
#' @param acomp_zero_impute Character ("CZM", "GBM","SQ","BL") or NULL; indicating weather to perform replacement of 0 abundance values with an estimate of the probability that the zero is not 0 (implemented only for avg_type = "acomp")
#' @param group Variable name in \code{\link[phyloseq]{sample_data}}) which defines sample groups for averaging (default is NULL)
#' @param drop_group_zero Logical; indicating weather OTUs with zero abundance withing a group of samples should be removed
#' @param progress Name of the progress bar to use ("none" or "text"; see \code{\link[plyr]{create_progress_bar}})
#' @param ... Additional arguments may be passed to \code{\link[zCompositions]{cmultRepl}}
#' @details
#' Typical OTU abundance tables in metagenomic analysis usually has different
#' sampling effort for different samples (which is an artifact of the sequencing procedure).
#' The total number of reads is meaningless and distance between OTU compositions is
#' on the relative scale (e.g., OTUs with 1 and 2 reads in one sample are so far as OTUs
#' with 10 and 20 reads in the other samples). Therefore such OTU tables represents
#' closed compositions and requires a special treatment within Aitchison geometry framework.
#'
#' Zero OTU abundance could be due to the insufficient number of reads. However,
#' it is possible to replace the zero counts with an expected value.
#' Bayesian-multiplicative (BM) replacement of count zeros is implemented in
#' \code{\link[zCompositions]{cmultRepl}} function of zCompositions package.
#' Sevral methods are supported: geometric Bayesian multiplicative (zero_impute = "GBM"),
#' count zero multiplicative (zero_impute = "CZM", default), Bayes-Laplace BM (zero_impute = "BL"),
#' or square root BM (zero_impute = "SQ").
#' In case of structural zeroes in OTU abundance table (e.g., absence of OTU
#' within a group assumes that it is not observed due to some biological pattern
#' and is not caused by a detection limit) "drop_group_zero" argument may be set
#' to "TRUE" to avoid zero replacement.
#'
#' @return phyloseq object with OTU relative abundance averaged over samples (all together or within a group).
#' @export
#' @seealso \code{\link[zCompositions]{cmultRepl}}, \code{\link[compositions]{acomp}}
#' @references
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V and Egozcue JJ (2017) Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 8:2224. doi: 10.3389/fmicb.2017.02224
#' Martin-Fernandez JA, Barcelo-Vidal C, Pawlowsky-Glahn V. (2003) Dealing With Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation. Mathematical Geology 35:3. doi: 10.1023/A:1023866030544
#'
#' @examples
#'
phyloseq_average <- function(physeq, avg_type = "acomp", acomp_zero_impute = NULL, 
    group = NULL, drop_group_zero = FALSE, verbose = TRUE, ...){

  # require(compositions)   # for Aitchison CoDa approach
  # require(zCompositions)  # for Bayesian-multiplicative replacement
  # require(plyr)

  ## Add warning for zero imputation and non-CoDa averaging
  if(!is.null(acomp_zero_impute) & avg_type != "acomp"){
    warning("Warning: imputations of zeros is implemented only for CoDa approach.\n")
    zero_impute <- NULL
  }

  ## Progress indicator
  if(verbose == TRUE){
    progress <- "text"
  } else {
    progress <- "none"
  }

  ## Average througth the all samples
  if(is.null(group)){
    res <- OTU_average(physeq, avg_type = avg_type, acomp_zero_impute = acomp_zero_impute, verbose = verbose)
  } ## End of single group

  ## Average by group
  if(!is.null(group)){

    ## Backup phylogenetic tree and remove it from the main object
    with_phy_tree <- phyloseq::phy_tree(physeq, errorIfNULL = F)
    if(!is.null(with_phy_tree)){
      tree <- phyloseq::phy_tree(physeq)
      physeq@phy_tree <- NULL
    }

    ## Split data by group
    ph_gr <- phyloseq_sep_variable(physeq, variable = group, drop_zeroes = drop_group_zero)

    ## Average OTU proportions within each group
    res <- plyr::llply(.data = ph_gr, .fun = OTU_average, avg_type = avg_type, acomp_zero_impute = acomp_zero_impute, verbose = verbose, .progress = progress)

    ## Give the group names to the averaged proportions
    for(i in 1:length(res)){
      phyloseq::sample_names(res[[i]]) <- names(res)[i]
    }

    ## Combine results
    # do.call(merge_phyloseq, res)   # doesn't work
    res_mrg <- res[[1]]
    for(j in 2:length(res)){
        res_mrg <- phyloseq::merge_phyloseq(res_mrg, res[[j]])
    }
    res <- res_mrg
    rm(res_mrg)

    ## Restore phylogenetic tree
    if(!is.null(with_phy_tree)){
      phyloseq::phy_tree(res) <- tree
    }
  } ## End of multiple groups

  ## Add averaging details as attributes to the resulting list
  attr(res, which = "Average_type") <- avg_type
  if(avg_type == "acomp"){ 
    attr(res, which = "Average_ZeroImputationMethod") <- acomp_zero_impute
  } else {
    attr(res, which = "Average_ZeroImputationMethod") <- NULL
  }
  attr(res, which = "Average_ByGroup") <- group
  attr(res, which = "Average_GroupZerosRemoved") <- drop_group_zero

  return(res)
}


## Function to average OTU relative abundances
OTU_average <- function(x, avg_type = "aldex", 
  acomp_zero_impute = NULL, aldex_samples = 128, aldex_denom = "all",
  result = "phyloseq", verbose = TRUE){
  # x = phyloseq object
  # avg_type = averaging type ("aldex" for ALDEx2-based averaging, "acomp" for Aitchison CoDa approach; "arithmetic" for simple arithmetic mean)
  # acomp_zero_impute = NULL or character indicating the method of zero imputation ("CZM" or "GBM")
  # aldex_samples = The number of Monte-Carlo Dirichlet instances to generate
  # aldex_denom = which features to use as the denominator for the geometric mean calculation ("zero", "all", "iqlr", "lvha")
  # result = resulting object ("phyloseq" or "matrix")
  # verbose = logical; shows warnings

  ## Remove sample metadata
  if(!is.null(phyloseq::sample_data(x, errorIfNULL = FALSE))){
    x@sam_data <- NULL
  }

  ## Extract OTU abundance table
  otus <- as.data.frame(phyloseq::otu_table(x))

  ## Show some warnings
  if(verbose == TRUE){
    # How many zeros are in the table? If more than 15% - rise a warning
    nz <- sum(otus == 0)
    if(nz > nrow(otus)*ncol(otus)*0.15){
      warning("Warning: there are more than 15% of zeroes in OTU table. Consider some additional data filtering.\n")
    }
    
    ## How many OTUs are with zero total abundance?
    oz <- phyloseq::taxa_sums(x) == 0
    if(any(oz)){
      warning("Warning: there are ", sum(oz), " OTUs with zero total abundance.\n")
    }
  }

  ## Transpose OTU abundance table (samples must be ROWS from this step!)
  if(phyloseq::taxa_are_rows(x) == TRUE){
    otus <- t(otus)
  }

  ## Acomp-based averaging
  if(avg_type == "acomp"){

      ## Replace 0 values with an estimate of the probability that the zero is not 0
      if(!is.null(acomp_zero_impute) & any(otus == 0)){
        otus <- try(
          zCompositions::cmultRepl(otus, label=0, method=acomp_zero_impute, output="prop", suppress.print = TRUE, ...)  # output="counts"  [zCompositions]
          )
        # Methods:
        #   CZM = multiplicative simple replacement
        #   GBM = Geometric Bayesian multiplicative
  
        if(class(otus) %in% "try-error"){
          stop("Error in multiplicative zero replacement, try to use other methods (e.g., 'zero_impute = CZM').\n")
        }
      }
  
      ## Transform the data using the the centred log-ratio (Aitchison compositions)
      otucomp <- compositions::acomp(otus)            # [compositions]
  
      ## Average proportions
      # TO DO: add possibilty to specify a robust estimator ('robust = TRUE')
      otuavg <- try( suppressWarnings( compositions::mean.acomp(otucomp, robust = F) ), silent = T)
      if(class(otuavg) %in% "try-error"){ stop("Error: method acomp failed.\n") }
      
      otuavg <- as.matrix(otuavg)    # it will be transposed here
      colnames(otuavg) <- "Average"  # rename average proporion column
  }

  ## ALDEx2-based averaging (thanks to Thom Quinn for the advice on these steps)
  if(avg_type == "aldex"){
    ## ALDEx2 requires OTUs as rows and samples as columns
    otus <- t(otus)

    ## Generate Monte Carlo samples of the Dirichlet distribution for each sample.
    ## Convert each instance using the centred log-ratio transform
    ald <- try( ALDEx2::aldex.clr(reads = otus, conds = rep("a", ncol(otus)), 
               mc.samples = aldex_samples, denom = aldex_denom, verbose = verbose, useMC = FALSE) )

    if(class(ald) %in% "try-error"){ stop("Error in ALDEx2::aldex.clr. Try to use another 'denom'.\n") }

    ## Extract the Monte Carlo Dirochlet instances
    mc <- ALDEx2::getMonteCarloInstances(ald)
    k <- ALDEx2::numMCInstances(ald)            # number of samples

    logratio <- 0   # placeholder

    ## Extract each Monte Carlo instance
    ## And add i-th log-ratio transformation to cumulative sum
    ## Based on propr code by Thom Quinn
    ## https://github.com/tpq/propr/blob/90a96e6e24397dcd50e89b25af00edc0fcf56197/R/aldex2propr.R#L70
    for(i in 1:k){
      mci_lr <- t(sapply(mc, function(x) x[, i]))
      logratio <- logratio + mci_lr
    }

    ## Average log-ratios
    logratio <- logratio / k

    ## Functions
    closure <- function(x){ x/sum(x) }               # divide proportions by sums
    geomean <- function(x){ prod(x)^(1/length(x)) }  # geometric mean

    ## Get the portions back
    ## ALDEx2 uses binary logarithm (log2) instead of natural logarithm, so for the back-transformation of log-ratios use 2^x instead of exp(x)
    mci_means <- t(apply(logratio, MARGIN = 1, FUN = function(x){ closure(2^x) }))

    ## Average OTU proportions across the samples
    otuavg <- apply(mci_means, MARGIN = 2, FUN = geomean)
    otuavg <- as.matrix(otuavg)
    colnames(otuavg) <- "Average"
  }

  ## Simple arithmetic averaging (in case if there are a lot of zeros and CoDa gives strange results)
  if(avg_type == "arithmetic"){

      ## Convert counts to proportions
      if(any(rowSums(otus) > 1)){  ## check that values are not proportions
         k <- .Machine$double.eps
         tmp <- pmax(k, apply(otus, MARGIN = 1, sum, na.rm = TRUE))
         otus <- sweep(otus, MARGIN = 1, tmp, "/")
      }

      ## Arithmetic averaging of proportions
      otuavg <- colMeans(otus, na.rm = TRUE)
      otuavg <- as.matrix(otuavg)
      colnames(otuavg) <- "Average"  # rename average proporion column
  }

  ## Back-transpose OTU abundances if neccesary
  # if(taxa_are_rows(x) == FALSE){
  #   otuavg <- t(otuavg)
  # }

  ## Replace original counts with the average relative abundance
  if(result == "phyloseq"){
    otu_table(x) <- phyloseq::otu_table(otuavg, taxa_are_rows = TRUE)
    return(x)
  }

  ## Just return matrix with averaged OTU proportions
  if(result == "matrix"){
    return(otuavg)
  }

}
