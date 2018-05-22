
## Prepare frequency counts for breakaway package
prepare_breakaway <- function(OTUs, remove_singletons = FALSE){

  ## If input is phyloseq or otu_table - extract OTU abundances
  if(any(c("phyloseq", "otu_table") %in% class(OTUs))){

    ## Check the orientation of the OTU table
    trows <- phyloseq::taxa_are_rows(OTUs)

    ## Extract OTUs
    if("phyloseq" %in% class(OTUs)){
      OTUs <- as.data.frame(phyloseq::otu_table(OTUs))
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
  res <- plyr::alply(.data = as.matrix(OTUs), .margins = 2, .fun = extract_non_zero)
  names(res) <- as.character(attr(res, "split_labels")[,1])


  ## Function to prepare data in breakaway format
    # index = an index variable
    # frequency = number of taxa that were observed with this frequency
  freq_count <- function(x, remove_singletons = FALSE){

    ## Remove singleton counts
    if(remove_singletons == TRUE){
      singls <- x == 1
      if(any(singls)){ x <- x[-which(singls)] }
    }

    rez <- table(x)
    rez <- as.data.frame(rez)
    colnames(rez) <- c("index", "frequency")
    return(rez)
  }

  ## Preprare frequency count table
  fct <- llply(.data = res, .fun = freq_count, remove_singletons = remove_singletons)

  return(fct)
}
