
#' @title Compare two phyloseq objects.
#'
#' @param phy1 Phyloseq object
#' @param phy2 Phyloseq object
#' @param cols Character vector with column names for phy1 & phy2 in the resulting table.
#'
#' @return Data frame with number and percentage of OTUs.
#' @export
#'
#' @examples
#' data(esophagus)
#'
#' # Remove taxa that have less than 5 reads
#' eso_trim <- prune_taxa(taxa_sums(esophagus) >= 5, esophagus)
#'
#' phyloseq_compare(esophagus, eso_trim, cols = c("Esophagus", "Trimmed esophagus"))
#'
phyloseq_compare <- function(phy1, phy2, cols = c("Before", "After")){

  # number of reads per OTU
  t1 <- taxa_sums(phy1)
  t2 <- taxa_sums(phy2)

  # percentage of reads
  if(sum(t1) > sum(t2)){ pr1 <- 100; pr2 <- (sum(t2)*100)/(sum(t1)) }
  if(sum(t1) < sum(t2)){ pr2 <- 100; pr1 <- (sum(t1)*100)/(sum(t2)) }

  # percentage of OTUs
  if(ntaxa(phy1) > ntaxa(phy2)){ po1 <- 100; po2 <- (ntaxa(phy2)*100)/(ntaxa(phy1)) }
  if(ntaxa(phy1) < ntaxa(phy2)){ po2 <- 100; po1 <- (ntaxa(phy1)*100)/(ntaxa(phy2)) }

  res <- rbind(
    data.frame(V0 = "Number of OTUs", V1 = ntaxa(phy1), V2 = ntaxa(phy2)),
    data.frame(V0 = "Percentage of OTUs", V1 = po1, V2 = po2),
    data.frame(V0 = "Total number of reads", V1 = sum(t1), V2 = sum(t2)),
    data.frame(V0 = "Percentage of reads", V1 = pr1, V2 = pr2),
    data.frame(V0 = "Average number of reads per OTU", V1 = mean(t1), V2 = mean(t2)),
    stringsAsFactors = F)

  colnames(res) <- c("Parameter", cols)
  return(res)
}
