
## Convert UC file (clustering output from USEARCH and VSEARCH) to Sumaclust OTU observation map
uc_to_sumaclust <- function(x, output_file = NULL){

  require(plyr)

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
