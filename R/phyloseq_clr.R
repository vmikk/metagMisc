
#' Centered log-ratio transformations for phyloseq objects
#'
#' Apply centered log-ratio (CLR) or modified CLR (mCLR) transformation to the
#' `otu_table` of a `phyloseq` object. Counts can be optionally shifted by a
#' user-specified pseudocount before converting to relative abundances and
#' transforming. Zeros are handled using a CLR implementation that ignores
#' zeros in the log and centers by the mean of the non-zero logs. The mCLR
#' variant applies an additional positive shift to non-zero CLR values to
#' ensure strict positivity, following Yoon et al. (2019).
#'
#' @param physeq A `phyloseq` object whose `otu_table` contains counts or
#'   abundances to be transformed.
#' @param method Character string selecting the transform: either
#'   "CLR" (default) for standard centered log-ratio with zero handling, or
#'   "mCLR" for the modified CLR that shifts non-zero CLRs to be strictly positive.
#' @param pseudocount Numeric scalar pseudocount to add to the `otu_table`
#'   prior to conversion to relative abundances. Set to `NULL` to disable.
#'   Default is 0.5.
#'
#' @return A `phyloseq` object whose `otu_table` has been transformed by the
#'   requested method.
#'
#' @details
#' The CLR implementation is adapted from the approach used in
#' `SpiecEasi` (Kurtz et al., 2015), computing `log(x)` only for entries above a
#' small tolerance and centering by the mean of the non-zero log values.
#'
#' The mCLR implementation follows Yoon et al. (2019) as used in the `SPRING`
#' package: CLR is computed on non-zero entries and then shifted by a positive
#' constant so that all non-zero values are strictly positive. In this
#' implementation, defaults are used (`eps = NULL`, `atleast = 1`).
#'
#' Prior to transformation, counts are converted to relative abundances per
#' sample using `phyloseq::transform_sample_counts`. If `pseudocount` is not
#' `NULL`, it is added to the counts before this conversion.
#'
#' @references
#' 
#' Kurtz ZD, Müller CL, Miraldi ER, Littman DR, Blaser MJ, et al. (2015) 
#' Sparse and Compositionally Robust Inference of Microbial Ecological Networks //
#' PLOS Computational Biology 11(5): e1004226. DOI:10.1371/journal.pcbi.1004226  
#' 
#' Yoon G, Gaynanova I, Müller CL (2019) 
#' Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data // 
#' Front. Genet. 10:516. DOI:10.3389/fgene.2019.00516  
#'
#' @seealso [phyloseq::transform_sample_counts()], [phyloseq::otu_table()]
#'
#'
#' @importFrom phyloseq otu_table transform_sample_counts
#' @export
#' 
phyloseq_clr <- function(physeq, method = "CLR", pseudocount = 0.5){

  ## CLR function with special handling of zeros
  # https://github.com/zdk123/SpiecEasi/blob/e4d7c0af49eae1ead40a89c1f569ee589de94776/R/normalization.R#L53
  clr_func <- function(x, base = exp(1), tol=.Machine$double.eps) {
    nzero <- (x >= tol)
    LOG <- log(ifelse(nzero, x, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }

  ## mCLR (Yoon et al., 2019)
  # https://github.com/GraceYoon/SPRING/blob/8049b2ae39adced84883d1ba3636c17a4de7180b/R/helpers.R#L61
  mclr_func <- function(x, base = exp(1), eps = NULL, atleast = 1, tol=.Machine$double.eps){
    # x = vector with values
    # eps = epsilon in eq (2) in "Yoon, Gaynanova, M\"{u}ller (2019), Frontiers in Genetics".
    #        positive shifts to all non-zero compositions. eps = absolute value of minimum of log ratio counts plus C
    # atleast = Constant C which ensures all nonzero values to be strictly positive (default is 1)

    ## CLR
    nzero <- (x >= tol)
    LOG <- ifelse(nzero, log(x, base), 0.0)                   # log only nonzero values, zeros stay as zeros
    clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0) # CLR

    ## Shift non-zero CLRs
    if(is.null(eps)){
      if(atleast < 0){
        warning("`atleast` should be positive. The functions uses default value 1 instead.\n")
        atleast = 1
      }
      if( min(clrdat) < 0 ){ # find the smallest negative value and add 1 to shift all data larger than zero.
        positivecst <- abs(min(clrdat)) + atleast # default `atleast` = 1
      } else {
        positivecst <- 0
      }
      # positive shift
      ADDpos <- ifelse(nzero, clrdat + positivecst, 0.0) ## make all non-zero values strictly positive.
      return(ADDpos)
    } else if(eps == 0){
      ## no shift. CLR transform applied to non-zero proportions only. without pseudo count.
      return(clrdat)
    } else if(eps > 0){
      ## use user-defined eps for additional positive shift.
      ADDpos <- ifelse(nzero, clrdat + eps, 0.0)
      return(ADDpos)
    } else {
      stop("Check your eps value for additional positive shift. Otherwise, leave it as NULL.\n")
   }
  }

  ## Add pseudocount if reqiured
  if(!is.null(pseudocount)){
    otu_table(physeq) <- otu_table(physeq) + pseudocount
  }

  ## Convert to relative abundances
  physeq <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

  ## CLR-transform OTU table
  if(method == "CLR"){ 
    otu_table(physeq) <- phyloseq::transform_sample_counts(otu_table(physeq), clr_func)
  }

  if(method == "mCLR"){ 
    otu_table(physeq) <- phyloseq::transform_sample_counts(otu_table(physeq), mclr_func)
  }

  return(physeq)
}
