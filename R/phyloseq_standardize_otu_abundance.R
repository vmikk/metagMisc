
#' @title Standardize OTU abundance table
#' @description Standardize phyloseq OTU table with with methods from \code{\link[vegan]{decostand}} from vegan package.
#' @param physeq A phyloseq-class object
#' @param method Standardization method
#' @param ... Additional parameters may be passed to vegan \code{\link[vegan]{decostand}} function
#' @return phyloseq object with standardized OTU table.
#' @seealso \code{\link[vegan]{decostand}}, \code{\link{phyloseq_transform_css}}, \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{physeq_transform_anderson_log}}
#' @details For implemented methods - see \code{\link[vegan]{decostand}}.
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
  rhea_depth = min( phyloseq::sample_sums(physeq) ), rhea_round = TRUE,
  ...){

  ## Method implemente in vegan
  vegan_methods <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")

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


  ## Replace old otu_table with the new one
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)

  return(physeq)
}
