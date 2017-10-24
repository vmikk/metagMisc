
#' @title Compare two phyloseq objects.
#'
#' @param phy1 A phyloseq-class object
#' @param phy2 A phyloseq-class object
#' @param cols Character vector with column names for phy1 & phy2 in the resulting table
#' @param more_stats Logical; if TRUE, some additional OTU abundance statistics will be calculated
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
phyloseq_compare <- function(phy1, phy2, cols = c("Before", "After"), more_stats = F){

  # number of reads per OTU
  t1 <- taxa_sums(phy1)
  t2 <- taxa_sums(phy2)

  # percentage of reads
  if(sum(t1) > sum(t2)){ pr1 <- 100; pr2 <- (sum(t2)*100)/(sum(t1)) }
  if(sum(t1) < sum(t2)){ pr2 <- 100; pr1 <- (sum(t1)*100)/(sum(t2)) }
  if(sum(t1) == sum(t2)){ pr1 <- 100; pr2 <- 100 }

  # percentage of OTUs
  if(ntaxa(phy1) > ntaxa(phy2)){ po1 <- 100; po2 <- (ntaxa(phy2)*100)/(ntaxa(phy1)) }
  if(ntaxa(phy1) < ntaxa(phy2)){ po2 <- 100; po1 <- (ntaxa(phy1)*100)/(ntaxa(phy2)) }
  if(ntaxa(phy1) == ntaxa(phy2)){ po1 <- 100; po2 <- 100 }

  ## Prepare resulting table
  res <- rbind(
    data.frame(V0 = "Number of samples", V1 = nsamples(phy1), V2 = nsamples(phy2)),
    data.frame(V0 = "Number of OTUs", V1 = ntaxa(phy1), V2 = ntaxa(phy2)),
    data.frame(V0 = "Percentage of OTUs", V1 = po1, V2 = po2),
    data.frame(V0 = "Total number of reads", V1 = sum(t1), V2 = sum(t2)),
    data.frame(V0 = "Percentage of reads", V1 = pr1, V2 = pr2),
    data.frame(V0 = "Average number of reads per OTU", V1 = mean(t1), V2 = mean(t2)),
    stringsAsFactors = F)

  colnames(res) <- c("Parameter", cols)

  ## Add additional abundance statistics
  if(more_stats == TRUE){

    ## Coefficient of quartile variation function
    cqv <- function(x){
      q1 <- quantile(x, probs = 0.25)
      q3 <- quantile(x, probs = 0.75)
      res <- (q3-q1)/(q3+q1)
      return(res)
    }

    ## Additional statistics
    adds <- rbind(
      data.frame(V0 = "Median number of reads per OTU", V1 = median(t1), V2 = median(t2)),
      data.frame(V0 = "Min total OTU abundance", V1 = min(t1), V2 = min(t2)),
      data.frame(V0 = "Q1 of total OTU abundance", V1 = quantile(t1, probs = 0.25), V2 = quantile(t2, probs = 0.25)),
      data.frame(V0 = "Q3 of total OTU abundance", V1 = quantile(t1, probs = 0.75), V2 = quantile(t2, probs = 0.75)),
      data.frame(V0 = "Max total OTU abundance", V1 = max(t1), V2 = max(t2)),
      data.frame(V0 = "Coefficient of quartile variation in OTU abundance", V1 = cqv(t1), V2 = cqv(t2)),
      stringsAsFactors = F)

    colnames(adds) <- c("Parameter", cols)

    ## Add it to the main table
    res <- rbind(res, adds, stringsAsFactors = F)
    rownames(res) <- NULL
  }

  return(res)
}
