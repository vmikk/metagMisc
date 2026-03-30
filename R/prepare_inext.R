
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

    # require(plyr)

    ## If input is phyloseq or otu_table - extract OTU abundances
    if(any(c("phyloseq", "otu_table") %in% class(OTUs))){

      ## Check the orientation of the OTU table
      trows <- phyloseq::taxa_are_rows(OTUs)

      ## Extract OTUs
      if("phyloseq" %in% class(OTUs)){
        OTUs <- phyloseq_otu_to_df(OTUs)
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

    OTUs <- dfRowName(x = OTUs, name = "OTU")
    setDT(OTUs)
    setcolorder(x = OTUs, neworder = "OTU")

    ## Data validation

    # Check data for empty samples
    samp_sums <- colSums(OTUs[, -1], na.rm = TRUE)
    if(any(samp_sums == 0)){
      to_rm <- names(samp_sums[ samp_sums == 0 ])
      warning("Empty samples (n=", length(to_rm), ") were removed from the data (samples with zero total abundance).\n")
      OTUs <- OTUs[, !..to_rm]
      rm(to_rm)
    }

    # Check data for negative entries
    has_negative <- OTUs[, 
      any(unlist(lapply(.SD, function(x) any(x < 0, na.rm = TRUE)))), 
      .SDcols = -1]
    if(has_negative){
      stop("There are negative values in the abundance data.\n")
    }


    # Remove zero and missing abundances
    trim_table <- function(x){
      x <- x[ Abundance > 0 ]
      if(any(is.na(x$Abundance))){
        warning("There are NA values in the abundance data, they will be excluded.\n")
        x <- x[ !is.na(Abundance) ]
      }
      return(x)
    }

    ## Convert to long format
    n_cll <- nrow(OTUs) * (ncol(OTUs) - 1)
    if(n_cll < 50000000){
      ## Reshape data in one pass
      OTUL <- melt(OTUs, id.vars = "OTU", variable.name = "SampleID", value.name = "Abundance")
      OTUL <- trim_table(OTUL)
    } else {
      ## Reshape data in chunks
      n_chunks <- data.table::fcase(
        n_cll <  9e7,                3L,
        n_cll >= 9e7 & n_cll < 5e8,  6L,
        n_cll >= 5e8 & n_cll < 5e9,  8L,
        n_cll >= 5e9 & n_cll < 5e10, 9L,
        n_cll >= 5e10,               12L)

      chunks <- chunk(x = colnames(OTUs)[-1], n = n_chunks)

      OTUL <- list()
      for(i in 1:n_chunks){
        cls <- c("OTU", chunks[[i]])
        OTUL[[i]] <- melt(OTUs[, ..cls], id.vars = "OTU", variable.name = "SampleID", value.name = "Abundance")
        OTUL[[i]] <- trim_table(OTUL[[i]])
        rm(cls)
      }

      OTUL <- rbindlist(OTUL)
    }

    rm(OTUs); invisible(gc())
    setorder(x = OTUL, SampleID, -Abundance, OTU)

    # Check if OTU counts are integers
    if(any(OTUL$Abundance != as.integer(OTUL$Abundance), na.rm = TRUE)){
      warning("There are non-integer values in the data, results may be meaningless.\n")
    }

    ## Prepare a list with OTU abundance vectors
    res <- OTUL[, .(abund = list(as.integer(sort(Abundance, decreasing = TRUE)))), by = SampleID]
    res <- setNames(res$abund, as.character(res$SampleID))

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
      res.cor <- plyr::llply(.data = res, .fun = single_est)
      return(res.cor)
    }
}
