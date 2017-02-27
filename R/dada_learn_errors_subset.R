#' DADA2 learning the error rates based on a subset of data
#' @description Function to learn the error rates from a subset of the data (for the big data workflow)
#'
#' @param filts character vector with the names of pre-filtered fastq or fastq.gz files
#' @param n integer, how many samples to choose for learning (will be ignored if samps != NULL)
#' @param samps character vector, which samples to take for learning (overrides n-argument)
#' @param multitr logical, enables multithreading
#' @param verbose logical, shows some additional information
#'
#' @return Matrix with estimated error rates that can be used with dada-function.
#' @export
#'
#' @references \url{https://benjjneb.github.io/dada2/bigdata.html}
#'
#' @seealso \code{\link{dada}}
#' @examples
#' library(dada2)
#' # filtF and filtR are character vectors with the full path to fastq.gz-files
#' errsF <- learn_errors_subset(filtF, n = 30)
#' errsR <- learn_errors_subset(filtR, samps = attr(errsF, "selected_samples"))

dada_learn_errors_subset <- function(filts, n = 25, samps = NULL, multitr = TRUE, verbose = TRUE){
    require(dada2)
    # require(plyr)

    if(is.null(samps)){
        # Choose N samples at random to learn from
        filts.learn <- sample(filts, n)
    }
    if(!is.null(samps)){
        # Take the proposed samples for learning
        filts.learn <- samps
    }

    # Dereplicate selected samples
    drp.learn <- derepFastq(filts.learn)

    if(verbose == TRUE){
        # Count number of reads in the selected samples
        readn <- plyr::laply(.data = drp.learn, .fun = function(z){ sum(z$uniques) })
        cat("There are ", sum(readn), " reads in the ", length(filts.learn), " selected samples.\n")
    }

    # Estimate error rates
    dd.learn <- dada(drp.learn, err=NULL, selfConsist=TRUE, multithread=multitr)
    err <- dd.learn[[1]]$err_out

    # Record which samples were selected for learning
    attr(err, "selected_samples") <- filts.learn

    return(err)
}
