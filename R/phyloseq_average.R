
#' @title Average relative OTU abundances.
#' @description This function implements OTU abundance averaging following CoDa (Compositional Data Analysis) workflow.
#' @param physeq A phyloseq-class object
#' @param zero_impute Character ("CZM", "GBM","SQ","BL") or logical; indicating weather to perform replacement of 0 abundance values with an estimate of the probability that the zero is not 0
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
#' @examples
#'
phyloseq_average <- function(physeq, zero_impute = "CZM", group = NULL, drop_group_zero = FALSE, progress = "text", ...){

  require(compositions)   # for Aitchison CoDa approach
  require(zCompositions)  # for Bayesian-multiplicative replacement
  require(plyr)

  ## If zero imputation method is specified as logical, sustitute it with CZM or NULL
  if(zero_impute == TRUE){
    zero_impute <- "CZM"
  }
  if(zero_impute == FALSE){
    zero_impute <- NULL
  }

  ## Function to average OTU relative abundances
  single_group_avg <- function(x, zeroimp = FALSE, meth = "CZM"){
    # x = phyloseq object
    # zeroimp = logical; if TRUE, zeros will be imputed
    # meth = method of zero imputation ("CZM" or "GBM")

    ## Remove sample metadata
    if(!is.null(sample_data(x, errorIfNULL = FALSE))){
      x@sam_data <- NULL
    }

    ## Extract OTU abundance table
    otus <- as.data.frame(otu_table(x))

    # How many zeros are in the table? If more than 15% - rise a warning
    nz <- sum(otus == 0)
    if(nz > nrow(otus)*ncol(otus)*0.15){
      warning("Warning: there are more than 15% of zeroes in OTU table. Consider some additional data filtering.\n")
    }

    ## Transpose OTU abundance table (samples must be ROWS from this step!)
    if(taxa_are_rows(x) == TRUE){
      otus <- t(otus)
    }

    ## Replace 0 values with an estimate of the probability that the zero is not 0
    if(zeroimp == TRUE){
      otus <- try(
        cmultRepl(otus, label=0, method=meth, output="prop", suppress.print = TRUE, ...)  # output="counts"  [zCompositions]
        )
      # Methods:
      #   CZM = multiplicative simple replacement
      #   GBM = Geometric Bayesian multiplicative

      if(class(otus) %in% "try-error"){
        stop("Error in multiplicative zero replacement, try to use other methods (e.g., 'zero_impute = CZM').\n")
      }
    }

    ## Transform the data using the the centred log-ratio (Aitchison compositions)
    otucomp <- acomp(otus)       # [compositions]

    ## Average proportions
    # TO DO: add possibilty to specify a robust estimator ('robust = TRUE')
    otuavg <- mean(otucomp)        # [compositions]
    otuavg <- as.matrix(otuavg)    # it will be transposed here
    colnames(otuavg) <- "Average"  # rename average proporion column

    ## Back-transpose OTU abundances if neccesary
    # if(taxa_are_rows(x) == FALSE){
    #   otuavg <- t(otuavg)
    # }

    ## Replace original counts with the average relative abundance
    otu_table(x) <- otu_table(otuavg, taxa_are_rows = TRUE)

    return(x)
  } ## End of single_group_avg


  ## Average througth the all samples
  if(is.null(group)){

    ## Without zero imputation
    if(is.null(zero_impute)){
      res <- single_group_avg(physeq, zeroimp = FALSE)
    }

    ## With zero imputation
    if(!is.null(zero_impute)){
      res <- single_group_avg(physeq, zeroimp = TRUE, meth = zero_impute)
    }
  } ## End of single group

  ## Average by group
  if(!is.null(group)){

    ## Backup phylogenetic tree and remove it from the main object
    with_phy_tree <- phy_tree(physeq, errorIfNULL = F)
    if(!is.null(with_phy_tree)){
      tree <- phy_tree(physeq)
      physeq@phy_tree <- NULL
    }

    ## Split data by group
    ph_gr <- phyloseq_sep_variable(physeq, variable = group, drop_zeroes = drop_group_zero)

    ## Average OTU proportions within each group
    if(is.null(zero_impute)){
      res <- llply(.data = ph_gr, .fun = single_group_avg, zeroimp = FALSE, .progress = progress)
    } else {
      res <- llply(.data = ph_gr, .fun = single_group_avg, zeroimp = TRUE, meth = zero_impute, .progress = progress)
    }

    ## Give the group names to the averaged proportions
    for(i in 1:length(res)){
      sample_names(res[[i]]) <- names(res)[i]
    }

    ## Combine results
    # do.call(merge_phyloseq, res)   # doesn't work
    res_mrg <- res[[1]]
    for(j in 2:length(res)){
        res_mrg <- merge_phyloseq(res_mrg, res[[j]])
    }
    res <- res_mrg
    rm(res_mrg)

    ## Restore phylogenetic tree
    if(!is.null(with_phy_tree)){
      phy_tree(res) <- tree
    }
  } ## End of multiple groups

  return(res)
}
