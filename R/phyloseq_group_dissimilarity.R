
#' @title Pairwise dissimilarity boxplots
#'
#' @param physeq A phyloseq-class object
#' @param group Grouping variable name (contained in \code{\link{sample_data}})
#' @param between_groups Logical, estimate between-group dissimilarity
#' @param method Distance/dissimilarity method (as character string; see \code{\link{distanceMethodList}})
#' @param method_title Logical, add method name to the plot title
#' @param notch Logical, draw the notches at each side of the boxes
#' @param justDF  Logical, instead of returning a ggplot2-object, return the data.frame that was used to build the plot
#' @param ... Additional arguments may be passed to \code{\link{distance}} function from phyloseq package
#'
#' @details If the notches of two boxplots do not overlap this indicates that the two medians differ.
#'
#' @return ggplot2-object or data.frame (if justDF = TRUE).
#'
#' @importFrom plyr llply ldply aaply
#' @importFrom ggplot2 ggplot aes geom_boxplot xlab ylab ggtitle
#' @export
#'
#' @examples
#' ## Load and subset data
#' data(enterotype)
#' ent <- subset_samples(enterotype, Enterotype %in% c("1", "2"))
#'
#' ## Dissimilarity boxplots
#' phyloseq_group_dissimilarity(ent, group = "Enterotype")
#' phyloseq_group_dissimilarity(ent, group = "Enterotype", between_groups = F)
#' phyloseq_group_dissimilarity(ent, group = "Enterotype", method = "jaccard")
#'
phyloseq_group_dissimilarity <- function(physeq, group = NULL, between_groups = TRUE, method = "bray", method_title = FALSE, notch = TRUE, justDF = FALSE, ...){

  ## Check the input
  if(is.null(sample_data(physeq, errorIfNULL = T))){ stop("Error: Sample data is missing in the phyloseq-object.\n") }
  if(is.null(group)){ stop("Error: groupping variable should be specified.\n") }

  ## Extract samle meta-data
  mtd <- as(object = sample_data(physeq), Class = "data.frame")
  mtd$SampleNames <- sample_names(physeq)

  if(!group %in% colnames(mtd)){
    stop("Error: Grouping variable is missing from the sample data of phyloseq-object.\n")
  }

  ## Count number of groups
  tabb <- table(mtd[, group])
  if(any(is.na(mtd[, group]))){ stop("Error: there are NA values in the grouping variable.\n") }
  if(length(tabb) == 1){ cat("Warning: there is only one group of samples in the resulting list.\n") }

  ## If there are several groups
  if(length(tabb) > 1){

    ## Split by groups
    physeq_split <- phyloseq_sep_variable(physeq, variable = group, drop_zeroes = T)

    ## Estimate within-group pairwise dissimilarities
    dd <- llply(
      .data = physeq_split,
      .fun = function(z, ...){ phyloseq::distance(z, method = method, type = "samples", ...) },
      ...)

    ## Convert dist to data frame
    ddm <- ldply(.data = dd, .fun = function(z){ 
      data.frame(Dist = as.vector(z), Comparison = "within-group")
      }, .id = "Group")

    ## Estimate between-group dissimilarities
    if(between_groups == TRUE){

      ## Full pair-wise comparisons
      ddf <- phyloseq::distance(physeq, method = method, type = "samples", ...)

      ## Convert dist to pairwise dataframe
      ddl <- dist2list(ddf, tri=FALSE)     # !! take both symmetric comparisons !! (a-b & b-a)

      ## Remove self-comparisons
      ddl <- ddl[-with(ddl, which(col == row)), ]

      ## Add group names to sample names
      ddl$ColGroup <- mtd[match(x = ddl$col, table = mtd$SampleNames), group]
      ddl$RowGroup <- mtd[match(x = ddl$row, table = mtd$SampleNames), group]

      ## Add same or different groups
      # ddl$SameGroup <- "no"
      # ddl$SameGroup[which(ddl$ColGroup == ddl$RowGroup)] <- "yes"

      ## Remove same-group-comparisons
      ddl <- ddl[which(!ddl$ColGroup == ddl$RowGroup), ]

      ## Sort group names alphabetically for each comparison (to take symmetric comparisons into account)
      grr <- t(apply(ddl[,c("ColGroup", "RowGroup")], 1, sort))

      ## Create comparison name (merge names of two groups)
      ddl$Group <- aaply(.data = grr, .margins = 1, .fun = paste, collapse = "-")

      ## Add to the main table
      ddm <- rbind(ddm, data.frame(Group = ddl$Group, Dist = ddl$value, Comparison = "between-groups"))

    } # end of between_groups

    ## Prepare plot
    pp <- ggplot(data = ddm, aes(x = Group, y = Dist, fill = Group)) +
            geom_boxplot(size = 0.8, notch = notch)
  } # end of two-groups case


  ## If only one group is provided
  if(length(tabb) == 1){

    ## Estimate pairwise dissimilarities
    dd <- phyloseq::distance(physeq, method = method, type = "samples", ...)

    ## Convert dist to data frame
    ddm <- data.frame(Dist = as.vector(dd), Group = names(tabb)[1])

    ## Plot pairwise dissimilarities
    pp <- ggplot(ddm, aes(x = Group, y = Dist)) +
            geom_boxplot(size = 0.8, notch = notch)

  } # end of one group

  ## Rename axes
  pp <- pp + xlab(group) +
    ylab(paste("Pairwise dissimilarity (", method, ")", sep=""))

  ## Add method name to the plot
  if(method_title == TRUE){
    pp <- pp + ggtitle(method)
  }

  if(justDF == FALSE){ return(pp) }   # return plot
  if(justDF == TRUE) { return(ddm) }  # return data frame with pairwise distances
}
