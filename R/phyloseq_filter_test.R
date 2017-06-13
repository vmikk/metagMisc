
## Automatic threshold selection for unabundat OTU trimming from phyloseq object
## based on a multivariate standard error estimation
phyloseq_filter_test <- function(physeq, group, dist_type = "bray", filter_by_group = TRUE, nresamp = 1000, ...){
    require(plyr)
    require(vegan)

    ## Extract grouping variable
    grp <- as(object = sample_data(physeq), Class = "data.frame")[, group]

    ## Filtering thresholds
    trh_prev <- expand.grid(prev.trh = seq(0.01, 0.20, 0.02), abund.trh = c(50, 100, 150, 200, 300))  ## NULL abund.trh will be tested sepatately
    trh_tota <- data.frame(frac = seq(0.001, 0.025, 0.004))                        # fraction of total abundance
    trh_rela <- data.frame(frac = c(0.00001, 0.0001, 0.0005, 0.001, 0.005, 0.01))  # relative abundance

    ## Apply different filtering thresholds
    batch_filter <- function(phys){
    
        ## Try-error filters
        phyloseq_filter_prevalence_tr <- function(...){
            rez <- try( phyloseq_filter_prevalence(...))
            if("try-error" %in% class(rez)){ rez <- NA }
            return(rez)
        }

        phyloseq_filter_taxa_tot_fraction_tr <- function(...){
            rez <- try( phyloseq_filter_taxa_tot_fraction(...))
            if("try-error" %in% class(rez)){ rez <- NA }
            return(rez)
        }

        phyloseq_filter_taxa_rel_abund_tr <- function(...){
            rez <- try( phyloseq_filter_taxa_rel_abund(...))
            if("try-error" %in% class(rez)){ rez <- NA }
            return(rez)
        }

        ## Filter low-prevalence OTUs
        # phyloseq_filter_prevalence(phys, prev.trh = 0.05, abund.trh = NULL) 
        f_prev <- mlply(.data = trh_prev, .fun = phyloseq_filter_prevalence_tr, physeq = phys)  # with abund.trh
        names(f_prev) <- paste("prev_",
            paste(
                "P_", attr(f_prev, "split_labels")[,"prev.trh"], ",",
                "A_", attr(f_prev, "split_labels")[,"abund.trh"],
                sep = ""), sep="")

        f_prevN <- mlply(.data = data.frame(prev.trh = unique(trh_prev[,"prev.trh"])),       # without abund.trh
            .fun = phyloseq_filter_prevalence_tr, physeq = phys, abund.trh = NULL)
        names(f_prevN) <- paste("prev_",
            paste("P_", attr(f_prevN, "split_labels")[,"prev.trh"], ",", "A_NULL", sep = ""),
            sep="")

        ## Remove OTUs with abundance less then a certain fraction of total abundance
        ## phyloseq_filter_taxa_tot_fraction(phys, frac = 0.01)  
        # f_tota <- mlply(.data = trh_tota, .fun = phyloseq_filter_taxa_tot_fraction_tr, physeq = phys)
        # names(f_tota) <- paste("rela_", attr(f_tota, "split_labels")[, "frac"], sep="")

        ## Remove OTUs with small mean relative abundance
        # phyloseq_filter_taxa_rel_abund(phys, frac = 1e-4)
        f_rela <- mlply(.data = trh_rela, .fun = phyloseq_filter_taxa_rel_abund_tr, physeq = phys)
        names(f_rela) <- paste("rela_", attr(f_rela, "split_labels")[, "frac"], sep="")

        ## Concatenate filtered data sets into one list
        rez <- c(f_prev, f_prevN,
                # f_tota,
                f_rela)

        ## Remove failed results (NA from try-error)
        f_na <- is.na(rez)
        if(any(f_na)){
            rez <- rez[-which(f_na)]
            attr(x = rez, which = "failed") <- names(f_na)  # add failed combination names
        }

        return(rez)
    }

    ## Filter data for all groups together
    if(filter_by_group == FALSE){
        cat("..Data filtering\n")
        RES <- batch_filter(physeq)
    }

    ## Independently filter data within each group
    if(filter_by_group == TRUE){

        ## Remove phylogenetic tree to speed up filtering
        if(!is.null(phy_tree(physeq, errorIfNULL = F))){
            physeq_no_tree <- physeq
            physeq_no_tree@phy_tree <- NULL
        } else {
            physeq_no_tree <- physeq
        }

        ## Split data by grouping variable
        cat("..Splitting data by grouping variable\n")
        physeq_grp <- phyloseq_sep_variable(physeq_no_tree, variable = group, drop_zeroes = TRUE)

        ## Batch filter OTUs within each group
        cat("..Data filtering\n")
        filt_grp <- llply(.data = physeq_grp, .fun = batch_filter, .progress = "text")

        ## Find which filtering thresholds failed and remove them
        failed_trh <- llply(.data = filt_grp, .fun = function(z){ attr(z, which = "failed") })
        fff <- unique( do.call(c, failed_trh) )
        
        ## Inform user about failed result
        if(!is.null(fff)){
            cat("Warning: Some filtering combinations failed with error and were excluded from the further consideration.\n")
            cat("See '...$failed_filters' for more details.\n")
        
            ## Remove failed thresholds from ALL groups
            filt_grp <- llply(.data = filt_grp, .fun = function(z){ 
                to_rm <- names(z) %in% fff
                if(any(to_rm)){
                    z <- z[-which(to_rm)]
                }
                return(z)})
        } # end of NULL fff
        
        ## Merge groups by filtered combination
        cmbs <- unique(unlist(lapply(filt_grp, names)))                           # all combination names
        RES <- list()
        for(i in cmbs){
            ## Extract the same combination for each group
            tmp <- llply(.data = filt_grp, .fun = function(z, trh = i){ z[[ which(names(z) == trh) ]] })

            ## Merge sample groups
            mrg <- tmp[[1]]
            for(j in 2:length(tmp)){
                mrg <- merge_phyloseq(mrg, tmp[[j]])
            }

            RES[[i]] <- mrg
            rm(tmp, mrg, j)
        }
        names(RES) <- cmbs

        ## Add original phylogenetic tree (if present)
        if(!is.null(phy_tree(physeq, errorIfNULL = F))){
            RES <- llply(.data = RES, .fun = function(z){ phy_tree(z) <- phy_tree(physeq); return(z) })
        }


    } # end of filter_by_group


    ## Function to estimate multSE for phyloseq
    calc_SE <- function(phys = physeq, grp = grp, dst = dist_type, ...){
    
        ## Create OTU-by-site distance matrix
        D <- phyloseq::distance(phys, method = dst, type = "samples")   
    
        ## Estimate multivariate standard error
        mse <- multSE(D, group = grp, ...)

        ## Perform PERMANOVA and extract pseudo F-ratio, R-squared and sums of squares
        adon <- vegan::adonis(D ~ grp, permutations = 1)
        mse <- data.frame(mse,
            F = adon$aov.tab["grp", "F.Model"],
            R2 = adon$aov.tab["grp", "R2"],
            SS.grp = adon$aov.tab["grp", "SumsOfSqs"],
            SS.tot = adon$aov.tab["Total", "SumsOfSqs"])

        return(mse)
    }

    ## Estimate multSE for initial dataset
    cat("..Estimating multSE for initial dataset\n")
    se0 <- calc_SE(physeq, grp, dist_type, progress = "none", nresamp = nresamp)

    ## Estimate multSE for each filtered dataset
    cat("..Estimating multSE for filtered datasets:\n")
    seF <- ldply(
        .data = RES,
        .fun = function(z){ calc_SE(z, grp, dist_type, progress = "none", nresamp = nresamp) },
        .progress = "text",
        .id = "FilteringThreshold")

    ## Merge results
    se_res <- rbind(
        data.frame(FilteringThreshold = "none", se0),
        seF)

    ## Estimate total abundance (number of reads that passed filtering) and it to the resulting table
    se_res$TotAbund <- c(sum(sample_sums(physeq)),
                ldply(.data = RES, .fun = function(z){ data.frame(TotAbund = sum(sample_sums(z))) })$TotAbund
                )

    ## Prepare results
    results <- list()
    results$se <- se_res             # resulting table with multSE
    results$filtered_data <- RES     # list with filtered phyloseq objects

    ## Add failed combinations to the attributes
    if(exists("failed_trh")){
        results$failed_filters <- failed_trh    # list with filtering thresholds that failed
    }

    return(results)
}
