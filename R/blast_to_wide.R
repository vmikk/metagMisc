# TO DO - add no_hit option (if !is.null(seqs)  --> fill with NAs)

## Convert BLAST results to a wide table
## For each query, matches should be sorted in the desired order
# max_hits = number of targets to preserve in a table
# taxonomy 
#    = NA, no taxonomy should be added to the table
#    = NULL, use default tax ranks (Kingdom - Species)
#    = custom character vector with column names that contain taxonomy information
# seqs = query sequences (XStringSet object from Biostrings package)
# refs = target sequences (XStringSet object from Biostrings package)
blast_to_wide <- function(blst, max_hits = 10, taxonomy = NULL, seqs = NULL, refs = NULL, verbose = TRUE) {

  ## Add query length to the table
  if(!is.null(seqs)){
    if(any(!blst$QueryName %in% names(seqs))){
      cat("Warning: not all BLAST queries are in the `seqs`.\n")
    }
    if(verbose == TRUE){ cat("..Adding query length to the table\n") }
    blst <- tibble::add_column(blst, QueryLength = Biostrings::width(seqs[blst$QueryName]), .before = "QueryStart")
  }
  
  ## Add reference length to the table
  if(!is.null(refs)){
    if(any(!blst$AccID %in% names(refs))){
      cat("Warning: not all BLAST targets are in the `refs`.\n")
    }
    if(verbose == TRUE){ cat("..Adding reference length to the table\n") }
    blst <- tibble::add_column(blst, TargetLen = Biostrings::width(refs[blst$AccID]), .after = "QueryLength")
  }
  
  ## Estimate query coverage = (query-to - query-from + 1)/ query-len
  if(!is.null(seqs)){
    if(verbose == TRUE){ cat("..Estimating coverage\n") }
    est_coverage <- function(x,y,len){ round( (max(c(x,y)) - min(c(x,y)) + 1)/ len , 3) }
    est_coverage <- Vectorize(est_coverage)
    
    qq <- with(blst, est_coverage(x = QueryEnd, y = QueryStart, len = QueryLength))
    if(verbose == TRUE){ cat("..Adding coverage to the table\n") }
    blst <- tibble::add_column(blst, QueryCoverage = qq, .before = "Evalue")
    rm(qq)
  }
  
  ## Enumerate hits with data.table
  setDT(blst)
  if(verbose == TRUE){ cat("..Enumerating BLAST hits\n") }
  blst <- blst[, HitNum := 1:.N, by = QueryName]
  if(verbose == TRUE){ cat("..Selecting top N hits\n") }
  blst <- blst[ HitNum <= max_hits ]
  
  
  ## Measure variables for melting
  varz <- c("AccID", "SeqIdentity", "AlignLen", "MismatchN", "GapOpenings", "QueryStart", "QueryEnd", "TargetStart", "TargetEnd")
  if(!is.null(seqs)){ varz <- c(varz, "QueryLength") }
  if(!is.null(refs)){ varz <- c(varz, "TargetLen") }
  if(!is.null(seqs)){ varz <- c(varz, "QueryCoverage") }
  varz <- c(varz, "Evalue", "BitScore")
  

  ## Add taxonomy to the table?
  if(!any(is.na(taxonomy))){
  
    ## If not tax ranks provided - use the default ones
    if(is.null(taxonomy)){ 
      taxonomy <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    }
  
    ## Remove missing tax ranks
    if(verbose == TRUE){ cat("..Subsetting taxonomy\n") }
    taxonomy <- taxonomy[ which(taxonomy %in% colnames(blst)) ]
  
    varz <- c(varz, taxonomy)
  }
  
  ## Reshape data to long format
  if(verbose == TRUE){ cat("..Reshaping data to long format\n") }
  blst_long <- melt(data = blst,
    id.vars = c("QueryName", "HitNum"),
    measure.vars = varz,
    variable.name = "VAR", value.name = "VAL", factorsAsStrings = FALSE)
  
  ## Reshape data to wide format (1 row = OTU)
  if(verbose == TRUE){ cat("..Reshaping data to wide format\n") }
  blst_wide <- dcast(
    data = blst_long,
    formula = QueryName ~ HitNum + VAR,
    value.var = "VAL")

  rm(blst_long)

  ## Convert to numeric
  if(verbose == TRUE){ cat("..Converting numeric data\n") }

  to_numer <- c("SeqIdentity", "AlignLen", "MismatchN", "GapOpenings", "QueryStart", 
  	"QueryEnd", "TargetStart", "TargetEnd", "QueryCoverage", "Evalue", "BitScore")
  if(!is.null(seqs)){ to_numer <- c(to_numer, "QueryLength", "QueryCoverage") }
  if(!is.null(refs)){ to_numer <- c(to_numer, "TargetLen") }

  to_numer_ids <- grepl(pattern = paste(to_numer, collapse="|"), x = colnames(blst_wide))
  to_numer_names <- colnames(blst_wide)[to_numer_ids]
  
  if("data.table" %in% class(blst_wide)){
    blst_wide[, (to_numer_names) := lapply(.SD, as.numeric), .SDcols = to_numer_names]
  } else {
    blst_wide[ to_numer_ids ] <- sapply(blst_wide[to_numer_ids], as.numeric)
  }

  ## Add sequences
  if(!is.null(seqs)){
    if(verbose == TRUE){ cat("..Adding sequences to the table\n") }
    SQDT <- data.table(OTU = names(seqs), Seq = as.character(seqs) )
    blst_wide <- merge(
      x = blst_wide,
      y = SQDT,
      by.x = "QueryName", by.y = "OTU", all.x = TRUE)
  }
  
  

  ############ Reorder columns in the results
  
  ## Universal columns
  if(verbose == TRUE){ cat("..Reordering columns\n") }
  col_order <- "QueryName"
  if(!is.null(seqs)){ col_order <- c(col_order, "Seq") }

  ## Hit-specific columns
  clz <- c("AccID", "SeqIdentity", "AlignLen", "MismatchN", "GapOpenings", 
    "QueryStart", "QueryEnd", "TargetStart", "TargetEnd")
  if(!is.null(seqs)){ clz <- c(clz, "QueryLength") }
  if(!is.null(refs)){ clz <- c(clz, "TargetLen") }
  if(!is.null(seqs)){ clz <- c(clz, "QueryCoverage") }
  clz <- c(clz, "Evalue", "BitScore")

  if(!any(is.na(taxonomy)) && !is.null(taxonomy)){
    clz <- c(clz, taxonomy)
  }

  clz <- paste(
    rep(1:max_hits, each = length(clz)),
    clz,
    sep = "_")
  
  
  ## Reorder columns in data.table
  col_order <- c(col_order, clz)
  setcolorder(blst_wide, col_order)

  if(verbose == TRUE){ cat("..All done\n") }
  return(blst_wide)
}