
#' @title Compare two phyloseq objects.
#'
#' @param physeq A phyloseq-class object
#' @param ... Optionally more phyloseq objects
#' @param cols Character vector with column names for phy1 & phy2 in the resulting table
#' @param more_stats Logical; if TRUE, some additional OTU abundance statistics will be calculated
#' @param long Logical; if TRUE, results will be returned in a long format
#'
#' @return Data frame with number and percentage of OTUs.
#' @details
#' Optionally (if more_stats = TRUE), some additional OTU abundance statistics may be estimated (min, max, median, Q1, Q3, and CQV of OTU abundance).
#' Coefficient of quartile variation (a.k.a., Quartile coefficient of dispersion) is preferred to the 'classical' coefficient of variation for the analysis of samples from nonnormal distributions.
#' It is estimated by the following formula: (Q3-Q1)/(Q3+Q1), where Q1 is the population 25th percentile and Q3 is the population 75th percentile.
#' @export
#'
#' @examples
#' data(esophagus)
#'
#' # Remove taxa that have less than 5 reads
#' eso_trim <- prune_taxa(taxa_sums(esophagus) >= 5, esophagus)
#'
#' phyloseq_compare(esophagus, eso_trim, cols = c("Esophagus", "Trimmed esophagus"))
#' phyloseq_compare(esophagus, eso_trim, cols = c("Esophagus", "Trimmed esophagus"), more_stats = T)
#'
phyloseq_compare <- function(physeq, ..., cols = NULL, more_stats = FALSE, long = FALSE){

  require(plyr)
  require(reshape2)

  ## Merge phyloseq objects into a list
  if(!missing(...)){
    lls <- list(physeq, ...)
  } else {   # if only one phyloseq was provided
    lls <- list(physeq)
  }

  ## Add names to each phyloseq object from the input
  if(is.null(cols)){
    names(lls) <- paste("Phys", 1:length(lls), sep="")
  } else {
    names(lls) <- cols
  }

  ## Summary function
  smr <- function(x, more_stats = FALSE){
    # x = single phyloseq object

    ## Number of reads per OTU
    treads <- taxa_sums(x)

    ## Number of reads per sample
    sreads <- sample_sums(x)

    ## Prepare resulting table
    res <- rbind(
      data.frame(V0 = "Number of samples", V1 = nsamples(x)),
      data.frame(V0 = "Number of OTUs", V1 = ntaxa(x)),
      data.frame(V0 = "Total number of reads", V1 = sum(treads)),
      data.frame(V0 = "Average number of reads per OTU", V1 = mean(treads)),
      data.frame(V0 = "Average number of reads per sample", V1 = mean(sreads)),
      stringsAsFactors = F)

    colnames(res) <- c("Parameter", "Phy")

    ## Add additional abundance statistics
    if(more_stats == TRUE){

      ## Coefficient of quartile variation function
      cqv <- function(x){
        q1 <- quantile(x, probs = 0.25)
        q3 <- quantile(x, probs = 0.75)
        rez <- (q3-q1)/(q3+q1)
        return(rez)
      }

      ## Additional statistics
      adds <- rbind(
        ## OTU-wise stats
        data.frame(V0 = "Median number of reads per OTU", V1 = median(treads)),
        data.frame(V0 = "Min total OTU abundance", V1 = min(treads)),
        data.frame(V0 = "Q1 of total OTU abundance", V1 = quantile(treads, probs = 0.25)),
        data.frame(V0 = "Q3 of total OTU abundance", V1 = quantile(treads, probs = 0.75)),
        data.frame(V0 = "Max total OTU abundance", V1 = max(treads)),
        data.frame(V0 = "Coefficient of quartile variation in OTU abundance", V1 = cqv(treads)),

        ## Sample-wise stats
        data.frame(V0 = "Median number of reads per sample", V1 = median(sreads)),
        data.frame(V0 = "Min total sample abundance", V1 = min(sreads)),
        data.frame(V0 = "Q1 of total sample abundance", V1 = quantile(sreads, probs = 0.25)),
        data.frame(V0 = "Q3 of total sample abundance", V1 = quantile(sreads, probs = 0.75)),
        data.frame(V0 = "Max total sample abundance", V1 = max(sreads)),
        data.frame(V0 = "Coefficient of quartile variation in sample abundance", V1 = cqv(sreads)),

        stringsAsFactors = F)

      colnames(adds) <- c("Parameter", "Phy")

      ## Add it to the main table
      res <- rbind(res, adds, stringsAsFactors = F)
      rownames(res) <- NULL
    } ## End of additional stats

    return(res)
  } ## End of summary function

  ## Apply summary function for each phyloseq object
  RES <- ldply(.data = lls, .fun = smr, .id = "Phyloseq", more_stats = more_stats)

  ## Reshape data to a wide format
  if(long == FALSE){
  
    ## Reshape
    RES <- reshape2::dcast(data = RES, Parameter ~ Phyloseq, value.var = "Phy")

    ## If there are multiple phyloseq objects - estimate percentages
    if(length(lls) > 1){
      
      ## Function to estimate percentage relative to the maximum value
      perc <- function(z){ 100 * z / max(z) }

      pct <- apply(
        X = rbind(
          RES[which(RES$Parameter == "Total number of reads"), -1],
          RES[which(RES$Parameter == "Number of OTUs"), -1]),
        MARGIN = 1, FUN = perc)

      pct <- data.frame(
        Parameter = c("Percentage of reads", "Percentage of OTUs"),
        t(pct))

      RES <- rbind(RES, pct)
    } ## End of percentages
  
  } ## End of long

  return(RES)
}
