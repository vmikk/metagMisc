## TO DO - add examples

#' @title Parse UC files (clustering output from USEARCH, VSEARCH, SWARM)
#'
#' @param x File name (typically with .uc extension)
#' @param map_only Logical, return only mapping (correspondence of query and cluster)
#' @param splitSeqID Logical, split sequence IDs at the semicolon. If TRUE (default), only the part of the ID before the semicolon is retained. If set to FALSE, the entire sequence ID, including any extra info after the semicolon, is kept intact
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
#'
#' If `rm_dups = TRUE`, duplicated sequences are removed based solely on their sequence IDs, without considering any size annotations that may be present.  
#'
#' @return Data.frame or data.table with the following columns:
#' \itemize{
#'   \item{"recordType"}{ - Indicates the type of the record (Hit, Centroid, Cluster, No hit).}
#'   \item{"clusterNumber"}{ - The number of the cluster, starting from 0.}
#'   \item{"seqLengthOrClusterSize"}{ - Represents either sequence length (for S, N, and H records) or cluster size (for C records).}
#'   \item{"percentIdentity"}{ - For H records, this shows the percent identity with the target sequence.}
#'   \item{"strand"}{ - For H records, indicates the DNA strand (+ or -) for nucleotides.}
#'   \item{"unusedField1"}{ - A placeholder for a field that is not used, kept for backward compatibility.}
#'   \item{"unusedField2"}{ - Another unused field, also included for backward compatibility.}
#'   \item{"alignment"}{ - CIGAR string, representing the compressed alignment, or an '=' sign indicating 100% identity to the target sequence.}
#'   \item{"queryLabel"}{ - Label of the query sequence.}
#'   \item{"targetLabel"}{ - Label of the target sequence.}
#' }
#' 
#' With `map_only = TRUE`, the resulting table will have only two columns - "Query" and "OTU".
#' 
#' @export
#' @import data.table
#'
#' @references
#' http://www.drive5.com/usearch/manual/opt_uc.html
#' @examples
#' parse_uc("usearch_OTUs.uc", map_only = F)
#' parse_uc("usearch_OTUs.uc", map_only = T)
#'
parse_uc <- function(x, map_only = FALSE, splitSeqID = TRUE, rm_dups = TRUE){

    ## Read the file
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

    ## Assign column names
    setnames(x = ii,
      old = paste0("V", 1:10),
      new = c(
        "recordType", "clusterNumber", "seqLengthOrClusterSize", 
        "percentIdentity", "strand",
        "unusedField1", "unusedField2",
        "alignment", "queryLabel", "targetLabel"),
      skip_absent = TRUE)

    ## Convert similarity to numeric
    ii[ percentIdentity %in% "*" , percentIdentity := NA ]  # There is no identity for centroids
    ii[ , percentIdentity := as.numeric(percentIdentity) ]

    ## Return only mapping results
    ## Subset to Query and OTU names only
    if(map_only == TRUE){
      if(splitSeqID == TRUE){
        ii <- ii[, .(Query, OTU)]
      } else {
        ii <- ii[, .(queryLabel, targetLabel)]
        setnames(x = ii,
          old = c("queryLabel", "targetLabel"),
          new = c("Query", "OTU"))
      }
    }

  return(ii)
}
