## TO DO: add example

#' @title Parse SWARM files
#' @description This function reads the results of SWARM clustering (*.swarms file).
#' @param x Input file name (usually with "swarm" extension)
#'
#' @details
#' SWARM is a single-linkage clustering method with iterative growth process and the use of sequence abundance values to delineate OTUs.
#' SWARM properly delineates large OTUs (high recall), and can distinguish OTUs with as little as two differences between their centers (high precision).
#'
#' @return data frame with two columns: OTU name and sequences names that belongs to this OTU.
#' @export
#' @seealso
#' Mahé F, Rognes T, Quince C, de Vargas C, Dunthorn M. (2015) Swarm v2: highly-scalable and high-resolution amplicon clustering. PeerJ 3:e1420 doi: 10.7717/peerj.1420
#' https://github.com/torognes/swarm
#' https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline
#'
#' @examples
#' parse_swarms("tst_1f.swarms")
#'
#'
parse_swarms <- function(x, otu_delimiter = " "){
    require(plyr)

    ## Load file
    sw <- readLines(x)

    ## Split swarm for each OTU
    sw <- alply(.data = sw, .margins = 1, .fun = function(z){ strsplit(z, split = otu_delimiter)[[1]] })

    ## Split abundance and OTU name
    sw <- llply(.data = sw, .fun = function(z){ do.call(rbind, strsplit(z, split = ";"))[,1] })

    ## Name OTUs as the first member of swarm
    sw.names <- vector()
    for(i in 1:length(sw)){
        sw.names <- c( sw.names, sw[[i]][1] )
    }
    names(sw) <- sw.names

    ## Convert result to data.frame
    sw <- ldply(.data = sw, .fun = function(z){ data.frame(Query = z, stringsAsFactors = F) }, .id = "OTU")
    sw$OTU <- as.character(sw$OTU)

    return(sw)
}
