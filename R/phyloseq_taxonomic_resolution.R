
#' @title Estimate taxonomic resolution of data
#'
#' @param physeq A phyloseq-class object
#' @param add_counts Logical, add read and OTU counts on top of the bars (default, TRUE)
#' @param justDF Logical, just return the data.frame, do not plot (default, FALSE)
#'
#' @return ggplot-object (if justDF = FALSE) or a data.frame (if justDF = TRUE).
#' @details This function counts the number of OTUs and reads annotated 
#' at the lowest level of taxonomic classification and shows it on a barplot.
#' @export
#'
#' @seealso \code{\link{get_max_taxonomic_rank}}
#' @examples
#' library(phyloseq)
#' library(ggplot2)
#' data(GlobalPatterns)
#' 
#' taxres <- phyloseq_taxonomic_resolution(GlobalPatterns)
#' taxres
#' taxres$data_wide
#' 
phyloseq_taxonomic_resolution <- function(physeq, add_counts = TRUE, justDF = FALSE){

  ## Determine the lowest level of taxonomic classification
  sp_ranks <- get_max_taxonomic_rank(physeq)

  ## Add number of reads for each taxa
  sp_ranks$NumReads <- phyloseq::taxa_sums(physeq)

  ## Count number of taxa vs rank
  sp_ranks_tab <- table(sp_ranks$RankName)

  ## Count number of reads vs rank
  sp_ranks_count <- plyr::ddply(
      .data = sp_ranks,
      .variables = "RankName",
      .fun = function(z){ data.frame(NumReads = sum(z$NumReads)) })

  ## Match number of taxa to number of reads
  sp_ranks_tab <- sp_ranks_tab[ match(x = sp_ranks_count$RankName, table = names(sp_ranks_tab)) ]

  ## Merge tables
  sp_ranks_ok <- data.frame(sp_ranks_count, Count = as.vector(sp_ranks_tab))
  sp_ranks_ok$PercNumReads <- with(sp_ranks_ok, NumReads / sum(NumReads))
  sp_ranks_ok$PercCount <- with(sp_ranks_ok, Count / sum(Count))

  ## Plot data
  if(justDF == FALSE){

    ## Reshape data
    sp_ranks_long <- reshape2::melt(
        data = sp_ranks_ok,
        id.vars = c("RankName", "NumReads", "Count"),
        variable.name = "DataType",
        value.name = "Perc")

    ## Create column with counts
    sp_ranks_long$Counts <- sp_ranks_long$NumReads
    sp_ranks_long$Counts[which(sp_ranks_long$DataType == "PercCount")] <- sp_ranks_long$Count[which(sp_ranks_long$DataType == "PercCount")]
    sp_ranks_long$NumReads <- NULL
    sp_ranks_long$Count <- NULL

    ## Make the plot
    pp <- ggplot(data = sp_ranks_long, aes(x = RankName, y = Perc, fill = DataType)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_fill_manual(
          values = c("#F6A600", "#E25822"),  # orangeyellow + reddishorange
          labels = c("Reads", "OTUs"),
          name = "") +
      labs(x = "Lowest taxonomic rank", y = "Percentage of data")

    ## Add absolute number of reads or OTUs to the bars
    if(add_counts == TRUE){
      pp <- pp + geom_text(
                    aes(label = Counts),
                    hjust = 0.5, vjust = -0.5,
                    size = 3.5,
                    position = position_dodge(0.9))
    }

    ## Add data in wide format (for export)
    pp$data_wide <- sp_ranks_ok

    return(pp)
  }

  ## Return only summary table
  if(justDF == TRUE){
    return(sp_ranks_ok)
  }

}
