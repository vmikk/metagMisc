
## Compute OTU diversity for a particular sample size or coverage
phyloseq_effort_div <- function(physeq, base = "size", level = NULL, conf = 0.95, correct_singletons = FALSE){
  # conf=0.95 (could be NULL)
  # If base="size" and level=NULL,
  #   then this function computes the diversity estimates for the minimum sample size among all sites.
  # If base="coverage" and level=NULL,
  #   then this function computes the diversity estimates for the minimum sample coverage among all sites.

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(phyloseq::otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Estimate species diversity
  res <- iNEXT::estimateD(x, datatype = "abundance", base = base, level = level, conf = conf)
  colnames(res)[1] <- "SampleID"

  ## Add sample metadata
  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){
    mtd <- as(phyloseq::sample_data(physeq), "data.frame")
    mtd$SampleID <- rownames(mtd)
    res <- merge(res, mtd, by = "SampleID")
  }

  return(res)
}
# phyloseq_effort_div(esophagus, level = 3000)
# phyloseq_effort_div(esophagus, level = 3000, conf = NULL)



#' @title Vizualize OTU diversity and richness for a particular sequencing depth or sample coverage ranges
#' @description
#' This function helps to estimate OTU diversity and sample coverage for a particular sequencing depth or
#' to estimate the required sequencing depth and OTU diversity for a particular sample coverage
#' @param physeq A phyloseq-class object
#' @param range Numeric vector of coverage levels or sequence depths to test, e.g., seq(from = 100, to = 1000, by = 100)
#' @param range_type "coverage" or "seqdepth"
#' @param variable Groupping variable name (contained in \code{\link[phyloseq]{sample_data}}); default, NULL
#' @param yvar Which dependent variable to show ("coverage", "seqdepth", "Q1", "Q2", "Q3").
#'   For range_type = "seqdepth": "coverage", "Q0", "Q1", "Q2".
#'   For range_type = "coverage": "seqdepth", "Q0", "Q1", "Q2"
#' @param interval Which interval to show on the group-averaged data ("CI", "quartiles", or NULL)
#' @param show_plot Logical; if TRUE, show the resulting plot
#' @param justDF Logical; if TRUE, just return the resulting table
#' @param progress Show progress bar ("text" or "none")
#' @param ... Additional arguments may be passed to \code{\link{phyloseq_effort_div}}
#' @details
#' Diversity estimation based on the analytical approach for coverage-based rarefaction proposed by Chao and Jost (2012).
#' @return ggplot object (if show_plot = TRUE) or a data.frame (if justDF = TRUE)
#' @export
#'
#' @examples
#' data("esophagus")
#' # Show sequencing depth required to reach a particular sample coverage
#' phyloseq_effort_div_rangeplot(esophagus, range = seq(from = 0.8, to = 0.97, by = 0.02),
#'     range_type = "coverage", yvar = "seqdepth", interval = "quartiles")
#' 
#' # Show OTU richness (Q = 0) achieved by different sequencing efforts
#' phyloseq_effort_div_rangeplot(esophagus, range = seq(from = 100, to = 500, by = 20),
#'     range_type = "seqdepth", yvar = "Q0", interval = "quartiles")
#' 
#' 
#' data("GlobalPatterns")
#' # Subset data
#' GP <- subset_samples(GlobalPatterns, SampleType %in% c("Soil", "Feces", "Skin"))
#' GP <- phyloseq_filter_top_taxa(GP, perc = 20)
#' 
#' # Show sequencing depth required to reach a particular sample coverage level for different sample groups
#' phyloseq_effort_div_rangeplot(GP, range = seq(from = 0.8, to = 0.97, by = 0.02),
#'   variable = "SampleType", range_type = "coverage", yvar = "seqdepth", interval = "quartiles")
#' 
#' # Show OTU richness (Q = 0) achieved by different sequencing efforts for different sample groups
#' phyloseq_effort_div_rangeplot(GP, range = seq(from = 100, to = 500, by = 20),
#'    variable = "SampleType", range_type = "seqdepth", yvar = "Q0")
#' 
phyloseq_effort_div_rangeplot <- function(physeq, range, range_type = "coverage",
                                          variable = NULL, yvar = "seqdepth", interval = "quartiles",
                                          show_plot = TRUE, justDF = FALSE, progress = "text", ...){

  # Q = Hill's q-value (default = 0 - OTU number; 1 = Shannon diversity) - in case of range_type = "seqdepth"

    ## TO DO - add automatic range sizes
  # range <- seq(from = 100, to = 500, by = 20)      # check seq depth distribution
  # range <- seq(from = 0.8, to = 0.98, by = 0.02)   # "coverage"

  ## Estimate required sequencing depth and OTU diversity for a particular sample coverage
  if(range_type == "coverage"){
    efforts <- plyr::mdply(
      .data = data.frame(level = range),
      .fun = function(...){ phyloseq_effort_div(physeq, base="coverage", conf = NULL, ...) },
      .progress = progress)

    ## Replace estimated coverage with the supplied coverage (they are different because of numerical issues)
    efforts$SC <- efforts$level
  }

  ## Estimate OTU diversity and sample coverage for a particular sequencing depth
  if(range_type == "seqdepth"){
    efforts <- plyr::mdply(
      .data = data.frame(level = range),
      .fun = function(...){ phyloseq_effort_div(physeq, base="size", conf = NULL, ...) },
      .progress = progress)
  }

  ## Return data (otherwise make a plot)
  if(justDF == TRUE){ return(efforts) }

  ## Which variables to take for plotting (colnames from iNEXT::estimateD results)
  if(yvar == "coverage"){ yy <- "SC"; ylab <- "Sample coverage" }
  if(yvar == "seqdepth"){ yy <- "m"; ylab <- "Sequencing depth (number of reads)" }
  if(yvar == "Q0"){ yy <- "q = 0"; ylab <- "Species richness, Q = 0" }
  if(yvar == "Q1"){ yy <- "q = 1"; ylab <- "Species diversity, Q = 1" }
  if(yvar == "Q2"){ yy <- "q = 2"; ylab <- "Species diversity, Q = 2" }

  if(range_type == "coverage"){ xx <- "SC"; xlab <- "Sample coverage" }
  if(range_type == "seqdepth"){ xx <- "m"; xlab <- "Sequencing depth (number of reads)" }

  if(interval == "CI"){ il <- "CI.lo"; iu <- "CI.up" }
  if(interval == "quartiles"){ il <- "Qart.lo"; iu <- "Qart.up" }

  ## Average results for plotting
  avgeff <- function(z, yy){
    z <- z[, yy]                                 # extract dependent variable
    cis <- ggplot2::mean_cl_boot(z, na.rm=TRUE)  # bootstap CI
    rez <- data.frame(
      Mean = cis[[1]],
      CI.lo = cis[[2]],
      CI.up = cis[[3]],
      Qart.lo = quantile(z, probs = 0.25),
      Qart.up = quantile(z, probs = 0.75))
    rownames(rez) <- NULL
    return(rez)
  }
  # e.g., avgeff(efforts, yy = "q = 1")

  ## If multiple groups are provided
  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE)) & !is.null(variable)){
    res_avg <- plyr::ddply(.data=efforts, .variables=c(variable, xx), .fun=avgeff, yy = yy)
    pp <- ggplot(data = res_avg, aes_string(x = xx, y = "Mean", color = variable, group = variable)) +
      geom_line() +
      geom_ribbon(aes_string(ymin = il, ymax = iu, color = NULL, fill = variable), alpha = 0.2)

  } else {

    ## Single group
    res_avg <- plyr::ddply(.data=efforts, .variables=xx, .fun=avgeff, yy = yy)
    pp <- ggplot(data = res_avg, aes_string(x = xx, y = "Mean")) +
      geom_line() +
      geom_ribbon(aes_string(ymin = il, ymax = iu), alpha = 0.2)
  }

  ## Add axes labels
  pp <- pp + labs(x = xlab, y = ylab)

  ## Add raw data to the results
  attr(pp, which = "Data") <- efforts
  attr(pp, which = "DataAverage") <- res_avg

  if(show_plot == TRUE){ print(pp) }
  invisible(pp)
}
