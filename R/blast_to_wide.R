# TO DO - add no_hit option (if !is.null(seqs)  --> fill with NAs)

#' @title Convert BLAST results to a wide table
#'
#' @param blst Data.table with blast hits (in long format, from m8 file)
#' @param max_hits Number of targets to preserve in a table
#' @param taxonomy \itemize{
#'   \item{"NA"}, no taxonomy should be added to the table
#'   \item{"NULL"}, use default tax ranks (Kingdom - Species)
#'   \item{"custom character vector"} with column names that contain taxonomy information
#' }
#' @param seqs Query sequences (XStringSet object from Biostrings package)
#' @param refs Target sequences (XStringSet object from Biostrings package)
#' @param add_coverage Logical, add query coverage to the output (default, TRUE)
#' @param remove_missing Logical, remove blast matches if sequences are missing in `seqs` and `refs` (default, FALSE)
#' @param verbose Logical, produce more screen output (default, TRUE)
#'
#' @details
#' For each query, matches should be sorted in the desired order (e.g., by BitScore or E-Value).
#'
#' @return data.table
#'
#' @export
#'
#' @examples
#'
blast_to_wide <- function(blst, max_hits = 10, 
  taxonomy = NULL, seqs = NULL, refs = NULL,
  add_coverage = TRUE, remove_missing = FALSE, verbose = TRUE) {

  ## Add query coverage
  if(add_coverage == TRUE){

    if(is.null(seqs)){
      cat("WARNING: query sequnces are not provided, coverage could not be estimated.\n")
    }
    if(is.null(refs)){
      cat("WARNING: reference sequnces are not provided, coverage could not be estimated.\n")
    }

    ## Add query length to the table
    if(!is.null(seqs)){
      
      query_in_seqs <- blst$QueryName %in% names(seqs)
      if(any(!query_in_seqs)){
        cat("WARNING: not all BLAST queries are in the `seqs`.\n")
        if(verbose == TRUE){ cat("..n = ", sum(!query_in_seqs), " sequences are missing\n") }

        if(remove_missing == TRUE){
          cat("..Removing BLAST queries that are missing in FASTA file ('seqs').\n")
          blst <- blst[ query_in_seqs, ]
        } else {
          stop("Please verify that all queries are present in the FASTA file ('seqs')\n")
        }

      }
      if(verbose == TRUE){ cat("..Adding query length to the table\n") }
      blst <- tibble::add_column(blst,
        QueryLength = Biostrings::width(seqs[blst$QueryName]),
        .before = "QueryStart")
    }
  
    ## Add reference length to the table
    if(!is.null(refs)){
      
      matches_in_refs <- blst$AccID %in% names(refs)
      if(any(!matches_in_refs)){
        cat("WARNING: not all BLAST targets are in the `refs`.\n")
        if(verbose == TRUE){ cat("..n = ", sum(!matches_in_refs), " records are missing,\n") }
        if(verbose == TRUE){ cat("..among them ", length(unique(blst$AccID[ !matches_in_refs ])), " unique accession IDs\n") }

        if(remove_missing == TRUE){
          cat("..Removing BLAST matches that are missing in the reference database ('refs' FASTA)\n")
          blst <- blst[ matches_in_refs, ]
        } else {
          stop("Please verify that all BLAST targets are present in the reference database ('refs' FASTA)\n")
        }

      }
      if(verbose == TRUE){ cat("..Adding reference length to the table\n") }
  
      if("AccID" %in% colnames(blst)){
        ## If target headers were split, use Accession ID from query sequences
        blst <- tibble::add_column(blst,
          TargetLen = Biostrings::width(refs[blst$AccID]),
          .after = "QueryLength")
      } else {
        ## Use original (full) target header
        blst <- tibble::add_column(blst,
          TargetLen = Biostrings::width(refs[blst$TargetName]),
          .after = "QueryLength")
      }
  
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
  
  } # end of add_coverage


  ## Enumerate hits with data.table
  setDT(blst)
  if(verbose == TRUE){ cat("..Enumerating BLAST hits\n") }
  blst <- blst[, HitNum := 1:.N, by = QueryName]
  if(verbose == TRUE){ cat("..Selecting top N hits\n") }
  blst <- blst[ HitNum <= max_hits ]


  ## Merge variables for melting
  varz <- c("SeqIdentity", "AlignLen", "MismatchN", "GapOpenings", "QueryStart", "QueryEnd", "TargetStart", "TargetEnd")
  if("AccID" %in% colnames(blst)){ varz <- c("AccID", varz) }  # if headers were split
  if(!is.null(seqs)){ varz <- c(varz, "QueryLength") }
  if(!is.null(refs)){ varz <- c(varz, "TargetLen") }
  if(!is.null(seqs)){ varz <- c(varz, "QueryCoverage") }
  varz <- c(varz, "Evalue", "BitScore")


  ## Add taxonomy to the table?
  if(!any(is.na(taxonomy))){

    ## If not tax ranks provided - use the default ones
    if(is.null(taxonomy)){
      if("AccID" %in% colnames(blst)){
        taxonomy <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      } else {
        taxonomy <- "TargetName"
      }
    }

    ## Remove missing tax ranks
    if(!"TargetName" %in% taxonomy){
        if(verbose == TRUE){ cat("..Subsetting taxonomy\n") }
        taxonomy <- taxonomy[ which(taxonomy %in% colnames(blst)) ]
    }
    varz <- c(varz, taxonomy)
  }

  ## Reshape data to wide format (1 row = OTU)
  if(verbose == TRUE){ cat("..Reshaping data to wide format\n") }

  blst_wide <- dcast(
    data = blst,
    formula = QueryName ~ HitNum,
    value.var = varz)


  ## Rename columns
  if(verbose == TRUE){ cat("..Renaming columns\n") }
  newnames <- data.table(OldName = colnames(blst_wide)[ ! colnames(blst_wide) %in% "QueryName" ])
  newnames[ , VarName := gsub(pattern = "_[0-9]+$", replacement = "", x = OldName) ]
  hitid <- sub(pattern = "_", replacement = "",
      x = regmatches(m = regexpr("_[0-9]+$", newnames$OldName), x = newnames$OldName))
  newnames[ , HitNum := as.integer(hitid) ]
  newnames[ , NewName := paste0(HitNum, "_", VarName) ]
  setnames(x = blst_wide, old = newnames$OldName, new = newnames$NewName)


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
  col_order <- col_order[ col_order %in% colnames(blst_wide) ]
  setcolorder(x = blst_wide, neworder = col_order)

  if(verbose == TRUE){ cat("..All done\n") }
  return(blst_wide)
}
