
## Estimate interpolated and extrapolated Hill numbers and build rarefaction curve
phyloseq_inext <- function(physeq, Q = 0, correct_singletons = FALSE, endpoint=NULL, knots = 40,
    multithread = FALSE, show_CI = TRUE, show_sample_labels = TRUE, show_plot = TRUE, justDF = FALSE, ...) {
  # physeq = phyloseq object
  # Q = Hill's q-value (default = 0 - OTU number; 1 = Shannon diversity)
  # correct_singletons = Logical; apply Good-Turing correction
  # endpoint = sample size for extrapolation (default = NULL = double reference sample size)
  # knots = number of equally-spaced sample sizes (x-axis) to analyze
  # multithread = Logical;
  # show_CI = Logical; show bootstap confidence interval on the plot
  # show_sample_labels = Logical; add sample labels to the plot
  # show_plot = Logical; show plot on screen
  # justDF = Logical; return table with rarefaction results and do not show the plot
  # ... = will be passed to iNEXT (e.g., conf = 0.95, nboot = 50)

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
    inext_res <- llply(
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
    res <- llply(.data = inext_res, .fun = function(z){ z$iNextEst })

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
    return(res)
  }

  ## Extract coordinates for sample labels
  samplabs <- ddply(.data = res, .variables = "SampleID", .fun = function(z){
    mid <- which.max(z$qD)
    rez <- data.frame(MaxQD = z[mid, "qD"], SampSize = z[mid, "m"])
    return(rez)
  })

  ## Split data to interpolated, observed & extrapolated parts
  resl <- dlply(.data = res, .variables = "method", .fun = function(z){ z })

  ## Prepare a plot
  pp <- ggplot(data = res, aes(x = m, y = qD, group = SampleID)) +  # color = color
    geom_line(data = resl$interpolated, linetype = "solid") +
    geom_line(data = resl$extrapolated, linetype = "dashed") +
    geom_point(data = resl$observed, size = 2)

  ## Show confinence interval
  if(show_CI == TRUE){
    pp <- pp +
      geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL, color = NULL), alpha = 0.2)   # fill = color
  }

  ## Show sample labels
  if(show_sample_labels == TRUE){
    pp <- pp +
      geom_text(data = samplabs, aes(x = SampSize, y = MaxQD, label = SampleID), size = 4, hjust = -0.5)     # color = color
  }

  ## Add axes labels
  pp <- pp + labs(x = "Sample Size", y = paste("Species diversity, q = ", Q, sep = ""))

  if(show_plot == TRUE){
    print(pp)
  }

  ## Add

  invisible(pp)
}
