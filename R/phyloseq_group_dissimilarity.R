
# Pairwise dissimilarity boxplots
phyloseq_pairwise_dissimilarity <- function(physeq, group = NULL, between_groups = TRUE, method = "bray", method_title = TRUE, notch = TRUE, ...){

  require(plyr)

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
  if(length(tabb) == 1){ cat("Warning: there is only one group of samples in the resulting list.\n") }
  if(length(tabb) > 2) { stop("Error: number of sample groups should not exceed 2.\n") }

  ## If two groups are provided
  if(length(tabb) == 2){

    ## Split by groups
    physeq_split <- phyloseq_sep_variable(physeq, variable = group, drop_zeroes = T)

    ## Estimate within-group pairwise dissimilarities
    dd <- llply(.data = physeq_split, .fun = function(z, ...){
        phyloseq::distance(z, method = method, type = "samples", ...) })

    ## Convert dist to data frame
    ddm <- ldply(.data = dd, .fun = function(z){ data.frame(Dist = as.vector(z)) }, .id = "Group")


    ## Estimate between-group dissimilarities
    if(between_groups == TRUE){

      ## Full pair-wise comparisons
      ddf <- phyloseq::distance(physeq, method = method, type = "samples", ...)

      ## Convert distance matrix to pairwised list
      dist2list <- function (dist, tri=TRUE) {
        if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }

        dat <- as.data.frame(as.matrix(dist))
        if (is.null(names(dat))) {
            rownames(dat) <- paste(1:nrow(dat))
        }
        value <- stack(dat)$values
        rnames <- rownames(dat)
        namecol <- expand.grid(rnames, rnames)
        colnames(namecol) <- c("col", "row")
        res <- data.frame(namecol, value)

        if(tri == TRUE){    # return only lower triangular part of dist
          res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
        }

        return(res)
      } # end of dist2list

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

      ## Create comparison name (merge names of two groups)
      ddl$Group <- paste(names(tabb), collapse = "-")

      ## Add to the main table
      ddm <- rbind(ddm, data.frame(Group = ddl$Group, Dist = ddl$value))

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


  ## Add method name to the plot
  if(method_title == TRUE){
    pp <- pp + ggtitle(method)
  }

  return(pp)
}

