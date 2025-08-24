
## CLR transform function
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


