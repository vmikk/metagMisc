
#' Summarize taxonomic richness, relative abundance, and occurrence by rank
#'
#' Computes summaries of OTU counts (alpha richness) and relative abundance 
#' aggregated at a specified taxonomic rank (e.g., Phylum), 
#' together with overall gamma richness across all samples.
#'
#' The function:
#' - counts number of OTUs per sample for each taxon with a specified rank via `phyloseq_ntaxa_by_tax()`;
#' - computes gamma richness on a combined-sample object via `phyloseq_combine_samples()`;
#' - aggregates counts by rank with `speedyseq::tax_glom()` and derives per-sample relative abundance (percent);
#' - estimates OTU occurrence per rank with `phyloseq_otu_occurrence()`;
#' - summarizes per-sample metrics (alpha richness and relative abundance) by their 
#'   median and median absolute deviation (MAD) for each taxon at the chosen rank.
#'
#' Output columns include the chosen rank column(s) and, for each taxon:
#' - `Alpha_NOTU_Median`, `Alpha_NOTU_MAD` - summaries of OTU counts per sample
#' - `Gamma_NOTU` - overall richness across all samples
#' - `Relabund_Median`, `Relaund_MAD` - summaries of relative abundance (in %) per sample
#' - `Occurrence` - fraction of samples where the taxon was present
#'
#' @param ps A `phyloseq` object
#' @param rnk Taxonomic rank to aggregate by (default, "Phylum")
#'
#' @return A `data.frame` with one row per taxon at the specified rank and
#'   summary columns as described above.
#' @seealso phyloseq_ntaxa_by_tax, phyloseq_combine_samples,
#'   phyloseq_standardize_otu_abundance, phyloseq_otu_occurrence,
#'   phyloseq_rename_with_tax, phyloseq_otu_to_df, dfRowName
#'
#' @importFrom plyr ddply join_all
#' @importFrom speedyseq tax_glom psmelt
#' @importFrom stats median mad
#' @importFrom phyloseq rank_names
#' @export
#'
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
