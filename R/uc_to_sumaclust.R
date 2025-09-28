
#' @title Convert UC file (clustering output from USEARCH and VSEARCH) to Sumaclust OTU observation map
#' @param x File name in UC-format (typically with .uc extension)
#' @param output_file Optional, character string with the name of the output file
#' @return This function invisibly returns a list, where each slot corresponds to a single OTU.
#' Each slot is a charater vector that contains sequence names that belongs to this OTU.
#' The first element of a vector is a name of the representative sequence of this OTU (OTU centoid)
#' which occurs twice in this vector (i.e., there are two identical sequence IDs in a vector
#' for OTUs that are represented by only one sequence).
#'
#' @importFrom plyr dlply laply
#' @export
#' @seealso \code{\link{parse_uc}}
#' @references https://git.metabarcoding.org/obitools/sumaclust/wikis/home
#' @examples
#' ## OTU clustering with Sumaclust
#' ## OTU observation map is written to 'sumaclust_map.txt'
#' # sumaclust -g -f -t 0.97 -O sumaclust_map.txt input_seqs.fasta
#'
#' ## OTU clustering with VSEARCH
#' ## Results are in USEARCH cluster format (UC) in 'vsearch_clusters.uc' file
#' # vsearch --cluster_fast input_seqs.fasta --id 0.97 --uc vsearch_clusters.uc
#'
#' ## Convert VSEARCH clustering results (in UC-format) to Sumaclust OTU observation map
#' # uc_to_sumaclust("vsearch_clusters.uc", output_file = "vsearch_to_sumaclust.txt")
#'
uc_to_sumaclust <- function(x, output_file = NULL){

  ## Read UC file
  uc <- parse_uc(x, map_only = TRUE)

  ## Prepare Sumaclust map (first value = centroid ID, then all sequence IDs)
  sumaclust_map <- function(z){
    rez <- c(z$OTU[1], z$Query)
    return(rez)
  }

  ## Batch
  res <- dlply(.data=uc, .variables="OTU", .fun=sumaclust_map, .progress="text")

  ## Write results to file
  if(!is.null(output_file)){
    resl <- laply(.data = res, .fun = function(z){ paste(z, collapse = "\t") })
    writeLines(resl, output_file)
  }

  invisible(res)
}
