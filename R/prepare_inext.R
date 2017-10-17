
#' @title Prepare a list of (corrected) abundance vectors (iNEXT input).
#' @description This function prepares input data for \code{\link[iNEXT]{iNEXT}} function from iNEXT package (Hsieh et al., 2016)
#' @param OTUs Data frame with OTU abundances (rows = species, cols = samples) or phyloseq/otu_table object
#' @param correct_singletons Logical, if TRUE singleton counts will be corrected with modified Good–Turing frequency formula (Chiu, Chao 2016)
#' @return List of abundance vectors (each element is a separate sample)
#' @export
#' @references
#' Hsieh TC, Ma KH, Chao A (2016) iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). Methods Ecol Evol, 7: 1451–1456. DOI 10.1111/2041-210X.12613
#'
#' Chiu C, Chao A (2016) Estimating and comparing microbial diversity in the presence of sequencing errors. PeerJ, 4: e1634. DOI 10.7717/peerj.1634
#'
#' @examples
#' # Load packages
#' if(!require(iNEXT)) { devtools::install_github("JohnsonHsieh/iNEXT"); library(iNEXT) }
#'
#' # Load data
#' data(esophagus)
#'
#' # Prepare input data for iNEXT
#' esophagus_inext <- prepare_inext(esophagus, correct_singletons = F)
#' esophagus_inext
#'
#' # Estimate rarefied and extrapolated number of species (Hill number with order q=0)
#' esophagus_q0 <- iNEXT(esophagus_inext)
#'
#' # Plot rarefaction curves
#' plot(esophagus_q0)
#'
prepare_inext <- function(OTUs, correct_singletons = T){

    require(plyr)

    ## If input is phyloseq or otu_table - extract OTU abundances
    if(any(c("phyloseq", "otu_table") %in% class(OTUs))){

      ## Check the orientation of the OTU table
      trows <- taxa_are_rows(OTUs)

      ## Extract OTUs
      if("phyloseq" %in% class(OTUs)){
        OTUs <- as.data.frame(otu_table(OTUs))
      }
      if("otu_table" %in% class(OTUs)){
        OTUs <- as.data.frame(OTUs)
      }

      ## Transpose OTU table (species should be arranged by rows)
      if(trows == FALSE){
        OTUs <- t(OTUs)
      }
    }

    ## Convert data to data frame
    if(!class(OTUs) %in% "data.frame"){
      OTUs <- as.data.frame(OTUs)
    }

    ## Data validation
    # Check data for empty samples
    samp_sums <- colSums(OTUs, na.rm = TRUE)
    if(any(samp_sums == 0)){
      warning("Empty samples were removed from the data (samples with zero total abundance).\n")
      OTUs <- OTUs[, -which(samp_sums == 0)]
    }

    # Check data for negative entries
    if(any(OTUs < 0, na.rm = TRUE)){
      stop("There are negative values in the abundance data.\n")
    }

    # Check if OTU counts are integers
    tmp_matr <- as.matrix(OTUs)
    if(!identical(all.equal(as.integer(tmp_matr), as.vector(tmp_matr)), TRUE)){
      warning("There are non-integer values in the data, results may be meaningless.\n")
    }
    rm(tmp_matr)


    ## Function to extract non-zero OTUs and sort OTU abundance
    extract_non_zero <- function(x) {
      rez <- sort(x[ which(x > 0) ], decreasing = T)
      rez <- as.integer(rez)
      return(rez)
    }

    ## Extract OTU abundances
    res <- alply(.data = as.matrix(OTUs), .margins = 2, .fun = extract_non_zero)
    names(res) <- as.character(attr(res, "split_labels")[,1])

    ## Just return vectors with OTU abundances if singleton correction is not required
    if(correct_singletons == FALSE){
      return(res)
    }

    # Correct the raw singleton counts with modified Good–Turing frequency formula
    # Chiu, Chao (2016) PeerJ 4:e1634
    if(correct_singletons == TRUE){

      ## This function is based on 'singleton.Est'
      single_est <- function(x){
        # x = vector of OTU abundance, e.g, x <- c(20, 10, 4, 4, 3, 2, 1, 1, 1, 1, 1)

        f2 <- sum(x == 2)  #doubleton
        f3 <- sum(x == 3)
        f4 <- sum(x == 4)
        f1 <- ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))

        singls <- which(x==1)
        if(length(singls) > 0) { x <- x[-singls] }
        if(f1 > 0){ x <- c(x, rep(1, round(f1))) }
        x <- x[x>0]
        x <- as.integer(x)
        return(x)
      }

      ## Apply correction for each abundance vector
      res.cor <- llply(.data = res, .fun = single_est)
      return(res.cor)
    }
}
