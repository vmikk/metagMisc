
# Taxonomy parsing function
parse_taxonomy_ufits <- function(char.vec){
    ## Convert ufits taxonomic ranks to QIIME-like style
    char.vec <- gsub(pattern = ":", replacement = "__", x = char.vec)
    char.vec <- gsub(pattern = ",", replacement = ";", x = char.vec)
    res <- parse_taxonomy_qiime(char.vec)
    return(res)
}


parse_ufits_tax <- function(x){
    # x = vector of character strings from ufits taxonomy

    ## Split OTUId to SequenceId (e.g., JQ976006) or MethodId (UTAX,SINTAX) + Taxonomy
    res <- strsplit(x, split = ";")
    res <- do.call(rbind, res)

    # Prepare list of taxonomic assignments
    res <- alply(.data = res[, 2], .margins = 1, .fun = parse_taxonomy_ufits)

    # Convert each vector to matrix
    res <- llply(.data = res, .fun = function(x){ t(as.matrix(x)) })

    # Prepare table
    res <- data.frame( do.call(rbind.fill.matrix, res), stringsAsFactors = F)

    # Sort columns
    tax.levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax.levels <- data.frame(Level = tax.levels, Rank = 1:length(tax.levels))

    mm <- match(x = colnames(res), table = tax.levels$Level)
    res <- res[, as.character(tax.levels$Level[sort(mm)]) ]

    return(res)
}

