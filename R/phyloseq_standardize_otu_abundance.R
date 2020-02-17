
#' @title Standardize OTU abundance table
#' @description Standardize phyloseq OTU table with with methods from \code{\link[vegan]{decostand}} from vegan package.
#' @param physeq A phyloseq-class object
#' @param method Standardization method
#' @param ... Additional parameters may be passed to vegan \code{\link[vegan]{decostand}} function
#' @return phyloseq object with standardized OTU table.
#' @seealso \code{\link[vegan]{decostand}}, \code{\link{phyloseq_transform_css}}, \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{physeq_transform_anderson_log}}
#' 
#' @details
#' Supported methods:
#' \itemize{
#' \item \strong{"total"} - convert data to relative abundances (divide by sample total);
#' \item \strong{"pa"} - convert OTU abundances to presence/absence scale (0/1);
#' \item \strong{"log"} - logarithmic transformation as suggested by Anderson et al. (2006): log (x) + 1 for x > 0, zeros are left as zeros, logarithm base = 2. Please note this is not log(x+1);
#' \item \strong{"hellinger"} - square root of method = "total" (Legendre & Gallagher 2001);
#' \item \strong{"max"} - divide by sample maximum;
#' \item \strong{"frequency"} - divide by sample total and multiply by the number of non-zero items, so that the average of non-zero entries is one;
#' \item \strong{"normalize"} - make sample sum of squares equal to one;
#' \item \strong{"range"} - standardize values into range 0 ... 1. If all values are constant, they will be transformed to 0;
#' \item \strong{"rank"} - replace abundance values by their increasing ranks leaving zeros unchanged. Average ranks are used for tied values;
#' \item \strong{"rrank"} - replace abundance values by relative ranks with maximum 1. Average ranks are used for tied values;
#' \item \strong{"standardize"} - scale OTU abundances within sample to zero mean and unit variance;
#' \item \strong{"wisconsin"} - Wisconsin double standardization where species are first standardized by maxima and then sites by site totals;
#' \item \strong{"chi.square"} - divide by sample sums and square root of OTU sums, and adjust for square root of matrix total (Legendre & Gallagher 2001). When used with the Euclidean distance, the distances should be similar to the Chi-square distance used in correspondence analysis.
#' }
#' 
#' For the implementation of "total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log" methods see \code{\link[vegan]{decostand}}.
#' 
#' @export
#'
#' @examples
#' # Load data
#' data("esophagus")
#'
#' # Total-sum scaling (TSS) normalization
#' phyloseq_standardize_otu_abundance(esophagus, method = "total")
#' # the same as
#' transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU) )
#' identical(
#'   phyloseq_standardize_otu_abundance(esophagus, method = "total"),
#'   transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU)) )
#'
#' # Presence-absence scaling (0/1)
#' phyloseq_standardize_otu_abundance(esophagus, method = "pa")
#'
#' # Logarithmic transformation as in Anderson et al., 2006
#' phyloseq_standardize_otu_abundance(esophagus, method = "log")
#'
#' # Hellinger standardization
#' phyloseq_standardize_otu_abundance(esophagus, method = "hellinger")
#'
phyloseq_standardize_otu_abundance <- function(physeq, method = "total",
  rhea_depth = NULL, rhea_round = TRUE, rhea_threshold = NULL,
  ...){

  ## Method implemente in vegan
  vegan_methods <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")
  other_methods <- c("rhea", "wisconsin")

  ## Check the orientation of the OTU table
  trows <- phyloseq::taxa_are_rows(physeq)
  if(trows == TRUE){ marg <- 2 } else { marg <- 1 }

  ## Extact OTU table
  comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")

  ## Standardize with vegan methods
  if(method %in% vegan_methods){

    ## Standardize community table with vegan
    comm_std <- vegan::decostand(comm, method, MARGIN = marg, ...)

  }

  ## Normalize counts via division to the sample size and then multiplication by the size of the smaller sample (or other seq depth)
  ## as implemented in "Rhea" package (Lagkouvardos et al., 2017, DOI:10.7717/peerj.2836; Reitmeier et al. 2020, DOI:10.21203/rs.2.21240/v1)
  if(method == "rhea"){

    ## If there is no user-provided fixed sample size - use minimum sampling size
    if(is.null(rhea_depth)){
      rhea_depth <- min( phyloseq::sample_sums(physeq) )
    }

    ## Convert data to relative abundances
    comm <- vegan::decostand(comm, method = "total", MARGIN = marg)

    ## Remove OTUs with low relative abundance
    if(!is.null(rhea_threshold)){

      comm[ comm < rhea_threshold ] <- 0

      ## Re-normalize by total abundance again (sample sums should be 1)
      comm <- vegan::decostand(comm, method = "total", MARGIN = marg)
    }

    ## Multiply by the desired samples size
    comm_std <- comm * rhea_depth

    ## Round abundance to the nearest integer below its current value
    # if(rhea_round == TRUE){ comm_std <- floor(comm_std) }

    ## Round abundance to the nearest integer while preserving the rounded sum
    if(rhea_round == TRUE){ 
      
      ## Function by josliber
      # https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
      smart_round <- function(x) {
        y <- floor(x)
        indices <- tail(order(x-y), round(sum(x)) - sum(y))
        y[indices] <- y[indices] + 1
        return(y)
      }

      comm_std <- as.data.frame(comm_std)
      comm_std <- apply(comm_std, MARGIN = marg, FUN = smart_round)
    }
  } # end of "rhea"

  ## Wisconsin Double standardization
  ## Species (MARGIN=2) are first standardized by maxima and then sites by site totals
  if(method == "wisconsin"){

    if(marg == 2){ mm <- c(1, 2) }
    if(marg == 1){ mm <- c(2, 1) }

    comm_std <- vegan::decostand(comm, "max", MARGIN = mm[1])
    comm_std <- vegan::decostand(comm_std, "tot", MARGIN = mm[2])

  }

  ## Replace old otu_table with the new one
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)

  return(physeq)
}
