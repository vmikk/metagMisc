## TO DO - add examples

#' @title Parse UC files (clustering output from USEARCH, VSEARCH, SWARM)
#'
#' @param x File name (typically with .uc extension)
#' @param map_only Logical, return only mapping (correspondence of query and cluster)
#' @param package Which package to use ("base" or "data.table")
#' @param rm_dups  Logical, remove duplicated entries (default, TRUE)
#'
#' @details USEARCH cluster format (UC) is a tab-separated text file.
#' Description of the UC file format (from USEARCH web-site: http://www.drive5.com/usearch/manual/opt_uc.html):
#' 1       Record type S, H, C or N (see table below).
#' 2       Cluster number (0-based).
#' 3       Sequence length (S, N and H) or cluster size (C).
#' 4       For H records, percent identity with target.
#' 5       For H records, the strand: + or - for nucleotides, . for proteins.
#' 6       Not used, parsers should ignore this field. Included for backwards compatibility.
#' 7       Not used, parsers should ignore this field. Included for backwards compatibility.
#' 8       Compressed alignment or the symbol '=' (equals sign). The = indicates that the query is 100% identical to the target sequence (field 10).
#' 9       Label of query sequence (always present).
#' 10      Label of target sequence (H records only).
#'
#' Record      Description
#' H       Hit. Represents a query-target alignment. For clustering, indicates the cluster assignment for the query. If ‑maxaccepts > 1, only there is only one H record giving the best hit. To get the other accepts, use another type of output file, or use the ‑uc_allhits option (requires version 6.0.217 or later).
#' S       Centroid (clustering only). There is one S record for each cluster, this gives the centroid (representative) sequence label in the 9th field. Redundant with the C record; provided for backwards compatibility.
#' C       Cluster record (clustering only). The 3rd field is set to the cluster size (number of sequences in the cluster) and the 9th field is set to the label of the centroid sequence.
#' N       No hit (for database search without clustering only). Indicates that no accepts were found. In the case of clustering, a query with no hits becomes the centroid of a new cluster and generates an S record instead of an N record.
#'
#' @return Data.frame or data.table.
#' @export
#' @import data.table
#'
#' @references
#' http://www.drive5.com/usearch/manual/opt_uc.html
#' @examples
#' parse_uc("usearch_OTUs.uc", map_only = F)
#' parse_uc("usearch_OTUs.uc", map_only = T)
#'
parse_uc <- function(x, map_only = F, package = "data.table", rm_dups = TRUE){

  ## Load data with built-in R commands  
  if(package %in% "base"){

    cat("Option to use the `base` package will be duplicated in future releas of metagMisc!\n")

    ## Read file
    ii <- read.delim(x, header = F, stringsAsFactors = F)

    ## Remove redundant S-records
    redund <- ii$V1 == "S"
    if(any(redund)){ ii <- ii[-which(redund), ] }

    ## Split Query name
    ii$Query <- do.call(rbind, strsplit(x = ii$V9, split = ";"))[,1]

    ## Split OTU name
    ii$OTU <- do.call(rbind, strsplit(x = ii$V10, split = ";"))[,1]

    ## OTU name = query name for centroids
    ii$OTU[which(ii$V1 == "C")] <- ii$Query[which(ii$V1 %in% c("S", "C"))]

    ## Check for duplicates
    tmp <- ii[ , c("Query", "OTU")]
    dups <- duplicated(tmp)
    if(any(dups)){
      cat("Warning: duplicated rows found!\n")

      ## Remove duplicated seqs
      if(rm_dups == TRUE){
        cat("..", sum(dups), " duplicates removed\n")
        ii <- ii[ - which(dups), ]
      }
    }

    if(map_only == TRUE){
        ii <- ii[, which(colnames(ii) %in% c("Query", "OTU"))]
    }
  } # end of `base` package


  ## Load data with `data.table` package
  if(package %in% "data.table"){
    ## Read file
    ii <- fread(file = x, header = FALSE, sep = "\t")

    ## Remove redundant S-records
    ii <- ii[ ! V1 %in% "S" ]

    ## Split Query name
    ii[, Query := tstrsplit(V9, ";", keep = 1) ]

    ## Split OTU name
    ii[, OTU := tstrsplit(V10, ";", keep = 1) ]

    ## OTU name = query name for centroids
    ii[ V1 %in% "C", OTU := Query ]

    ## Check for duplicates
    if(nrow(ii[, .(Query, OTU)]) != nrow(unique(ii[, .(Query, OTU)]))){
      cat("Warning: duplicated rows found!\n")

      ## Remove duplicated seqs
      if(rm_dups == TRUE){
        dups <- duplicated(ii[, .(Query, OTU)])
        cat("..", sum(dups), " duplicates removed\n")
        ii <- ii[ ! dups ]
      }
    }

    ## Subset to Query - OTU names only
    if(map_only == TRUE){
      ii <- ii[, .(Query, OTU)]
    }
  }

  return(ii)
}
