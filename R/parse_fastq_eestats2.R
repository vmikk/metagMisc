#' @title Load USEARCH fastq_eestats2 report
#' @description Load USEARCH fastq_eestats2 report showing how many reads will pass an expected error filter at different length thresholds.
#' @param file the name of the input file.
#' @param long a logical value indicating wheather the results should be returned in long format.
#'
#' @return data frame
#' @export
#'
#' @details 
#' Additional details on 'fastq_eestats2' function can be found on http://drive5.com/usearch/manual/cmd_fastq_eestats2.html
#' 
#' @examples
#' ## Run USEARCH
#' # usearch -fastq_eestats2 reads.fq -output eestats2.txt -length_cutoffs 200,300,10
#' parse_fastq_eestats2("eestats2.txt", long = F)
#' parse_fastq_eestats2("eestats2.txt", long = T)
#' 
#' library(ggplot2)
#' ee <- parse_fastq_eestats2("eestats2.txt", long = T)
#' ggplot(data = ee, aes(x = Length, y = NumberOfReads, group = MaxEE, color = as.factor(MaxEE))) + 
#'   geom_point(size = 2) + 
#'   geom_line() + 
#'   theme_classic()
#'
parse_fastq_eestats2 <- function(file, long = FALSE){

    # read.delim("R1.eestats2.txt", stringsAsFactors = FALSE, skip = 2)
    # data.table::fread(input = "R1.eestats2.txt")

    ## Load file
    inp <- readLines(file)

    ## Remove empty rows
    inp <- inp[-c(1:3,5)]

    ## Extract column names
    header <- strsplit(inp[1], split = "  +")[[1]]

    ## Prepare columns for % or reads
    header <- rep(header, times=c(1, rep(2, length(header)-1)))  # duplicate columns
    header[seq(3, length(header), by = 2)] <- paste(header[seq(3, length(header), by = 2)], ",%", sep="")

    ## Prepare main table
    res <- gsub(pattern = "\\(", replacement = " ", x = inp[2:length(inp)])
    res <- gsub(pattern = "%)", replacement = " ", x = res[2:length(res)])
    res <- gsub("^\\s+|\\s+$", "", res)   # trim leading and ending whitespaces
    res <- plyr::adply(.data = res, .margins = 1, .fun = function(z){ strsplit(z, split = " +")[[1]] }, .id = NULL)

    colnames(res) <- header
    res <- colwise(as.numeric)(res)

    ## Reshape data
    if(long == TRUE){

      ## Columns with percentage of reads
      perc.col <- grep(pattern = "%", x = colnames(res))
      
      ## Reshape only read counts
      res_subset <- res[,-perc.col]
      setDT(res_subset)
      mm <- melt(data = res_subset,
                 id.vars = header[1],
                 variable.name = "MaxEE", value.name = "NumberOfReads")
      setDF(mm)

      ## Rename
      mm$MaxEE <- gsub(pattern = "\\(No EE cutoff\\)", replacement = Inf, x = mm$MaxEE)
      mm$MaxEE <- gsub(pattern = "MaxEE ", replacement = "", x = mm$MaxEE)
      mm$MaxEE <- as.numeric(mm$MaxEE)

      res <- mm
    }

    return(res)
}

