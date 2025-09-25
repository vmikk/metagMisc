phyloseq_tax_summary <- function(ps, rnk = "Phylum"){

  ## Count number of OTUs by taxonomic rank (per sample)
  otus <- phyloseq_ntaxa_by_tax(ps, TaxRank = rnk, add_meta_data = FALSE)

  ## Count number of OTUs by taxonomic rank (overall, gamma richness)
  gammadiv <- phyloseq_combine_samples(ps)
  gmm <- phyloseq_ntaxa_by_tax(gammadiv, TaxRank = rnk, add_meta_data = FALSE)
  gmm$Sample <- NULL
  colnames(gmm)[which(colnames(gmm) %in% "N.OTU")] <- "Gamma_NOTU"

  ## Count number of reads by taxonomic rank
  glomed <- speedyseq::tax_glom(physeq = ps, taxrank = rnk, NArm = TRUE)
  readsr <- phyloseq_standardize_otu_abundance(glomed, method = "total")
  readsr <- speedyseq::psmelt(readsr)
  readsr$Abundance <- readsr$Abundance * 100

  ## Estimate OTU occurrence
  occ <- phyloseq_otu_occurrence(glomed)
  occ <- phyloseq_rename_with_tax(occ, taxrank = rnk)
  occ <- phyloseq_otu_to_df(occ)
  occ <- dfRowName(occ, name = rnk)

  ## Average by sample
  my_avg <- function(x){
    data.frame(
      # Mean = mean(x, na.rm=TRUE),
      Median = median(x, na.rm=TRUE),
      MAD = mad(x, na.rm=TRUE)
      # Min = min(x, na.rm=TRUE),
      # Max = max(x, na.rm=TRUE),
      # Q1 = quantile(x, probs = 0.25, na.rm = TRUE),
      # Q3 = quantile(x, probs = 0.75, na.rm = TRUE),
      # SD = sd(x, na.rm=TRUE),
      # N = length(x)
      )
  }

  otus_avg <- ddply(.data = otus,
    .variables = colnames(otus)[colnames(otus) %in% rank_names(ps)],
    .fun = function(z){ my_avg(z$N.OTU) })

  reads_avg <- ddply(.data = readsr,
    .variables = colnames(readsr)[colnames(readsr) %in% rank_names(ps)],
    .fun = function(z){ my_avg(z$Abundance) })

  grp_colz <- which(colnames(otus_avg) %in% rank_names(ps))
  colnames(otus_avg)[-grp_colz] <- paste("Alpha_NOTU_", colnames(otus_avg)[-grp_colz], sep = "")

  grp_colz <- which(colnames(reads_avg) %in% rank_names(ps))
  colnames(reads_avg)[-grp_colz] <- paste("Relabund_", colnames(reads_avg)[-grp_colz], sep = "")


  # res <- cbind(otus_avg, reads_avg[-grp_colz])
  res <- plyr::join_all(list(otus_avg, gmm, reads_avg, occ))
  res <- cbind(
    res[, rank_names(ps)],
    res[, which(!colnames(res) %in% rank_names(ps))])

  return(res)
}
