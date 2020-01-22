
## Replace missing taxonomy
phyloseq_taxonomy_imputation <- function(phys,
  unknown_taxon = "_unidentified", unknown_sp = " sp",
  make_unique = FALSE){

  ## If input is of class 'phyloseq'
  inp_class <- class(phys)
  if("phyloseq" %in% inp_class || "taxonomyTable" %in% inp_class){

    if(is.null(phyloseq::tax_table(phys, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }

    x <- as.data.frame(phyloseq::tax_table(phys), stringsAsFactors = F)
  } else {
    stop("Error: input data should be of class `phyloseq` or `taxonomyTable`")
  }

  ## Function to replace NAs with higher taxa names + unident string
  replace_col <- function(x, col_num = 2, unident = "unidentified"){
    nas <- is.na(x[, col_num])
    if(any(nas) == FALSE){
      return(x)              # no missing values
    } else {
      ## which values to replace
      to_repl <- which(nas)

      if(col_num == 1){
      	x[to_repl, col_num] <- unident
      } else {
        x[to_repl, col_num] <- paste(x[to_repl, col_num - 1], unident, sep = "")
      }
      return(x)
    }
  }

  clz <- colnames(x)
  sp_in_ranks <- clz %in% c("Species", "species", "sp")
  if(any(sp_in_ranks)){
    sp_id <- which(sp_in_ranks)
    non_sp_id <- (1:length(clz))[-sp_id]
  } else {
    sp_id <- NA
    non_sp_id <- 1:length(clz)
  }

  ## Replace higher taxa
  for(i in non_sp_id){
    x <- replace_col(x, col_num = i, unident = unknown_taxon)
  }

  ## Replace species names
  if(!is.na(sp_id)){
    x <- replace_col(x, col_num = sp_id, unident = unknown_sp)
  }

  ## Function to remove multiple unident strings (e.g., "_unidentified_unidentified")
  replace_unidents <- function(tt, strr = "_unidentified", spp = " sp"){
    ## tt = character vector
  
    ## Prepare regex for paterns
    # Multiple string occurrences
    multpatt <- paste("(", strr, ")(\\1+)", sep = "")
    unsp <- paste(strr, spp, sep = "")

    ## Replace "_unidetified_unidentified" with single occurrence
    # gsub(x = tt, pattern = "(_unidentified)(\\1+)", replacement = "_unidentified", perl = T)
    rez <- gsub(x = tt, pattern = multpatt, replacement = strr, perl = T)

    ## Replace "_unidetified sp" with " sp"
    rez <- gsub(x = rez, pattern = unsp, replacement = spp, perl = T)    

    return(rez)
  }

  ## Remove multiple unident strings from all tax columns
  x <- sapply(x, replace_unidents)

  ## Make species names unique
  if(make_unique == TRUE){
    x[, ncol(x)] <- base::make.unique(names = x[, ncol(x)], sep = ".")
  }

  ## Replace tax_table with the modified one
  phyloseq::tax_table(phys) <- phyloseq::tax_table(as.matrix(x))
  return(phys)
}


