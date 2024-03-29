
#' @title Parse taxonomy string from AMPtk (for single OTU)
#'
#' @param x Character vector of length 1 with taxonomic annotation of single OTU/species.
#' @param custom_ranks Named vector with customized prefixes for taxonomic ranks
#' @details This function will split
#' @return Named character vector with taxonomic ranks of OTU or species.
#' @export
#' @seealso \code{\link{parse_taxonomy_qiime}}
#' @references AMPtk: Amplicon ToolKit for NGS data (formally UFITS): https://github.com/nextgenusfs/amptk
#' @examples
#' parse_taxonomy_amptk("k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella,s:Mortierella parvispora")
#' parse_taxonomy_amptk("k:Fungi,p:Ascomycota,g:Chalara")
#'
#' # With customized ranks
#' parse_taxonomy_amptk(
#'  x = "do:Eukaryota,su:Amorphea,cd:Obazoa,gr:Opisthokonta,sg:Holomycota(Nucletmycea),ki:Fungi,sk:Dikarya,ph:Ascomycota,sd:Saccharomycotina,cl:Saccharomycetes,sc:Saccharomycetidae,or:Saccharomycetales,fa:Saccharomycetaceae,ge:Saccharomyces,sp:Saccharomyces cerevisiae,st:Saccharomyces cerevisiae var boulardii",
#'  custom_ranks = c(do = "Domain", su = "Supergroup", cd = "Clade", gr = "Group", sg = "Subgroup", ki = "Kingdom", sk = "Subkingdom", ph = "Phylum", sd = "Subdivision", cl = "Class", sc = "Subclass", or = "Order", fa = "Family", ge = "Genus", sp = "Species", st = "Strain"))
parse_taxonomy_amptk <- function(x, custom_ranks = NULL){

    ## Use default tax ranks
    if(is.null(custom_ranks)){
      Tranks <- c(
        d = "Domain",
        k = "Kingdom",
        p = "Phylum",
        c = "Class",
        o = "Order",
        f = "Family",
        g = "Genus",
        s = "Species")
    } else {
    ## Alternatively use the provided tax ranks
      Tranks <- custom_ranks
    }

    ### Old function
    ## Convert ufits taxonomic ranks to QIIME-like style
    # x <- gsub(pattern = ":", replacement = "__", x = x)
    # x <- gsub(pattern = ",", replacement = ";", x = x)
    # res <- phyloseq::parse_taxonomy_qiime(x)

    ### Updated function (without phyloseq dependency)
    
    ## Split string by tax ranks
    x <- strsplit(x, ",", TRUE)[[1]]
    
    ## Remove leading and trailing spaces
    x <- gsub(pattern = "^[[:space:]]{1,}", replacement = "", x = x)
    x <- gsub(pattern = "[[:space:]]{1,}$", replacement = "", x = x)


    if(length(x) > 0){

      if(x[1] == "No hit"){
        res <- ""
      } else {

        ## Extract tax rank prefixes
        ranks_prefs_found <- grep("^[[:alpha:]]{1,}:", x)  # match prefixes, e.g., "k:"
        ## If no prefixes found - return dummy ranks
        if(length(ranks_prefs_found) == 0L){
          warning("No taxonomic ranks prefixes were found (e.g., 'k:'). Dummy ranks will be added.\n")
          names(x) <- paste("Rank", 1:length(x), sep = "")
          res <- x
        } else {
        ## Otherwise, remove prefixes and add tax ranks as names
          res <- gsub(pattern = "^[[:alpha:]]{1,}:", replacement = "", x = x)
          rnks <- regmatches(m = regexpr("^[[:alpha:]]{1,}:", x), x = x)
          rnks <- gsub(pattern = ":", replacement = "", x = rnks)
          repranks <- Tranks[rnks]
          names(res)[ranks_prefs_found[!is.na(repranks)]] <- repranks[!is.na(repranks)]
        }
      }

    } else {
      warning("Empty taxonomy vector encountered.")
      res <- x
    }

  return(res)
}



#' @title Parse multiple taxonomy strings from AMPtk.
#'
#' @param x Vector of character strings from AMPtk taxonomy
#' @param withID Logical; set to TRUE (default) if charter strings contain sequence IDs
#' @param multithread Logical; run the function in parallel
#' @details
#' This function splits AMPtk (ex-UFITS) result of taxonomic annotation by taxonomic rank and combines it into one table.
#' The first part of taxonomy string is assumed to be an ID of the best hit (e.g., EF040844) or a taxonomic classificator (UTAX or SINTAX) followed by a semicolon.
#'
#' @return Data frame with OTUs or species as rows and their taxonomic ranks as columns.
#' @references AMPtk: Amplicon ToolKit for NGS data (formally UFITS): https://github.com/nextgenusfs/amptk
#' @seealso \code{\link{parse_taxonomy_amptk}} \code{\link{parse_taxonomy_qiime}}
#' @export
#' @examples
#' ## Taxonomy with sequence ID
#' tax <- c(
#' "UDB011192;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Russulaceae,g:Russula,s:Russula chloroides",
#' "SINTAX;k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Rhizophydiales",
#' "EF040844;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Piloderma",
#' "EU240039;k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella",
#' "UTAX;k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Amphinema,s:Amphinema byssoides",
#' "SINTAX;k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Spizellomycetales"
#' )
#'
#' parse_taxonomy_amptk_batch(tax)
#'
#' ## Taxonomy without sequence ID
#' tax2 <- c(
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Russulales,f:Russulaceae,g:Russula,s:Russula chloroides",
#' "k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Rhizophydiales",
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Piloderma",
#' "k:Fungi,p:Zygomycota,o:Mortierellales,f:Mortierellaceae,g:Mortierella",
#' "k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Atheliales,f:Atheliaceae,g:Amphinema,s:Amphinema byssoides",
#' "k:Fungi,p:Chytridiomycota,c:Chytridiomycetes,o:Spizellomycetales"
#' )
#'
#' parse_taxonomy_amptk_batch(tax2, withID = FALSE)
#'
parse_taxonomy_amptk_batch <- function(x, withID = TRUE, multithread = FALSE){
    
    # require(plyr)

  ###############
  ############### Cluster setup
  ###############

  ## Progress bar type for single-threaded plyr functions
  progr <- "text"
  parall <- FALSE

  ## Check if foreach and doParallel packages are available
  if(multithread){
    if(!requireNamespace("foreach", quietly = TRUE)){ stop("foreach package required for parallel plyr operation.\n") }
    if(!requireNamespace("doParallel", quietly = TRUE)){ stop("doParallel package is required.\n") }
  }

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

    ## plyr arguments for parallel execution
    progr <- "none"
    parall <- TRUE

    ## Disable load balancing
    paropts <- list(preschedule=TRUE)

    ## Start the cluster
    cl <- parallel::makeCluster(ncores)

    ## Register the parallel backend
    doParallel::registerDoParallel(cl)

    ## Load packages on cluster nodes
    parallel::clusterEvalQ(cl, library("phyloseq"))

    ## Send useful objects to the workers
    vars <- c("parse_taxonomy_amptk")
    parallel::clusterExport(cl=cl, varlist=vars, envir = environment())
  }


  ###############
  ############### Parsing
  ###############

  if(withID == TRUE){
    ## Split OTUId to SequenceId (e.g., JQ976006) or MethodId (UTAX,SINTAX) + Taxonomy
    res <- strsplit(x, split = ";")
    res <- do.call(rbind, res)

    # Prepare list of taxonomic assignments
    res <- plyr::alply(.data = res[, 2], .margins = 1, 
      .fun = parse_taxonomy_amptk, .progress = progr, .parallel = parall)
  } else {
    # Prepare list of taxonomic assignments
    res <- plyr::alply(.data = x, .margins = 1,
      .fun = parse_taxonomy_amptk, .progress = progr, .parallel = parall)
  }

  # Convert each vector to matrix
  res <- plyr::llply(.data = res, .fun = function(x){ t(as.matrix(x)) })

  # Prepare table
  res <- data.frame( do.call(rbind.fill.matrix, res), stringsAsFactors = F)

  # Sort columns
  if(ncol(res) > 1){
    tax.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax.levels <- data.frame(Level = tax.levels, Rank = 1:length(tax.levels))
  
    mm <- match(x = colnames(res), table = tax.levels$Level)
    res <- res[, as.character(tax.levels$Level[sort(mm)]) ]
  }

  return(res)
}
