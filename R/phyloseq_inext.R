
#' @title Estimate interpolated and extrapolated Hill numbers and sample coverage and construct rarefaction curve.
#'
#' @param physeq A phyloseq-class object
#' @param Q Diversity order of Hill number; 0 for species richness (default), 1 for Shannon diversity; 2 for Simpson diversity
#' @param curve_type Which data to show on the plot: "diversity" (default) or sample "coverage"
#' @param correct_singletons Logical; apply Good-Turing correction for singleon counts
#' @param endpoint Sample size for extrapolation (default = NULL, which corresponds to a double reference sample size for each sample)
#' @param knots Number of equally-spaced sample sizes (x-axis) to analyze (default, 40)
#' @param multithread Logical; if TRUE, attempts to run the function on multiple cores
#' @param show_CI Logical; show bootstap confidence interval on the plot
#' @param show_sample_labels Logical; add sample labels to the plot
#' @param show_plot Logical; show plot on screen
#' @param justDF Logical; return table with rarefaction results and do not show the plot
#' @param add_raw_data Logical; add attributes with iNEXT results to the output (see \code{\link[iNEXT]{iNEXT}})
#' @param ... Additional arguments may be passed to \code{\link[iNEXT]{iNEXT}} (e.g., conf = 0.95, nboot = 100)
#' @return Plot of class 'ggplot' or 'data.frame' (if 'justDF = TRUE')
#' @export
#' @seealso \code{\link[iNEXT]{iNEXT}}, \code{\link{prepare_inext}}
#' @examples
#'
phyloseq_inext <- function(physeq, Q = 0, curve_type = "diversity",
    correct_singletons = FALSE, endpoint=NULL, knots = 40,
    multithread = FALSE, show_CI = TRUE, show_sample_labels = TRUE,
    show_plot = TRUE, justDF = FALSE, add_raw_data = TRUE, ...) {

  ## Prepare a list of OTU abundance vectors
  x <- prepare_inext(
        as.data.frame(phyloseq::otu_table(physeq, taxa_are_rows = T)),
        correct_singletons = correct_singletons)

  ## Run rarefaction for all samples in a single thread
  if(multithread == FALSE){
    ## Estimate interpolated and extrapolated Hill numbers
    inext_res <- iNEXT::iNEXT(x, q = Q, datatype = "abundance", endpoint = endpoint, knots = knots, ...)

    ## Extract results (list with data.frames for each sample)
    res <- inext_res$iNextEst
  }

  ## Run rarefaction in parallel (separately for each sample)
  if(multithread == TRUE){

    ## Check if foreach and doParallel packages are available
    if(!requireNamespace("foreach", quietly = TRUE)){ stop("foreach package required for parallel plyr operation.\n") }
    if(!requireNamespace("doParallel", quietly = TRUE)){ stop("doParallel package is required.\n") }

    ## Check the platform type and specify the number of cores to use
    if(multithread && .Platform$OS.type == "unix"){
      ncores <- parallel::detectCores()
      if(is.numeric(multithread)){ ncores <- multithread }
      if(is.na(ncores) | is.null(ncores)){ ncores <- 1 }
    } else {
      ncores <- 1
      if(multithread && .Platform$OS.type=="windows"){
        warning("Multithreading has been DISABLED, as forking is not supported on .Platform$OS.type 'windows'.\n")
      }
    }

    ## Setup cluster
    if(ncores > 1){
      ## Disable load balancing
      paropts <- list(preschedule=TRUE)

      ## Start the cluster
      cl <- parallel::makeCluster(ncores)

      ## Register the parallel backend
      doParallel::registerDoParallel(cl)

      ## Load packages on cluster nodes
      # parallel::clusterEvalQ(cl, library("iNEXT"))

      ## Send useful objects to the workers
      vars <- c("x", "Q", "endpoint")
      parallel::clusterExport(cl=cl, varlist=vars, envir = environment())
    } else {
      multithread <- FALSE
    }

    ## Estimate interpolated and extrapolated Hill numbers
    inext_res <- plyr::llply(
        .data = x,
        .fun = function(z, ...){ iNEXT::iNEXT(x = z, datatype = "abundance", ...) },
        .parallel = multithread,
        q = Q, endpoint = endpoint, knots = knots, ...)   # additional iNEXT arguments may be passed here

    ## Stop the cluster
    if(ncores > 1){
      parallel::stopCluster(cl)
      rm(cl)
    }

    ## Extract results (list with data.frames for each sample)
    res <- plyr::llply(.data = inext_res, .fun = function(z){ z$iNextEst })

  } # end of multithread = TRUE


  ## Prepare single table for all samples
  res <- Map(cbind, res, SampleID = names(res))   # Add sample name to the table
  res <- do.call(rbind, res)                      # Merge list into a single data.frame
  rownames(res) <- NULL

  ## Add sample metadata
  if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){
    mtd <- as(phyloseq::sample_data(physeq), "data.frame")
    mtd$SampleID <- rownames(mtd)
    res <- merge(res, mtd, by = "SampleID")
  }

  ## Just return the table (otherwise, make a plot)
  if(justDF == TRUE){

    ## Add raw data to the results
    if(add_raw_data == TRUE){
      attr(res, which = "iNEXT") <- inext_res
    }

    return(res)
  } # end of justDF


  ## Extract coordinates for sample labels
  samplabs <- plyr::ddply(.data = res, .variables = "SampleID", .fun = function(z){
    mid <- which.max(z$qD)
    rez <- data.frame(SampSize = z[mid, "m"], MaxQD = z[mid, "qD"], MaxSC = z[mid, "SC"])
    return(rez)
  })

  ## Split data to interpolated, observed & extrapolated parts
  resl <- plyr::dlply(.data = res, .variables = "method", .fun = function(z){ z })

  ## Which variables to plot?
  if(curve_type == "diversity"){ YY <- "qD"; YYL <- "qD.LCL"; YYU <- "qD.UCL"; YYM <- "MaxQD"; ylab <- paste("Species diversity, q = ", Q, sep = "") }
  if(curve_type == "coverage") { YY <- "SC"; YYL <- "SC.LCL"; YYU <- "SC.UCL"; YYM <- "MaxSC"; ylab <- "Sample coverage" }

  ## Prepare a plot
  pp <- ggplot(data = res, aes_string(x = "m", y = YY, group = "SampleID")) +  # color = color
    geom_line(data = resl$interpolated, linetype = "solid") +
    geom_line(data = resl$extrapolated, linetype = "dashed") +
    geom_point(data = resl$observed, size = 2)

  ## Show confinence interval
  if(show_CI == TRUE){
    pp <- pp +
      geom_ribbon(aes_string(ymin = YYL, ymax = YYU, color = NULL), alpha = 0.2)   # fill = color
  }

  ## Show sample labels
  if(show_sample_labels == TRUE){
    pp <- pp +
      geom_text(data = samplabs, aes_string(x = "SampSize", y = YYM, label = "SampleID"), size = 4, hjust = -0.5)     # color = color
  }

  ## Add axes labels
  pp <- pp + labs(x = "Sample Size", y = ylab)

  if(show_plot == TRUE){
    print(pp)
  }

  ## Add raw data to the results
  if(add_raw_data == TRUE){
    attr(pp, which = "RarefTable") <- res
    attr(pp, which = "iNEXT") <- inext_res
  }

  invisible(pp)
}
