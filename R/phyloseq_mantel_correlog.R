
## Function to perform Mantel correlogram test
phyloseq_mantel_correlog <- function(ps,
  otu_dist = "bray",
  geo_dist = "geodesic", log_geo_dist = TRUE,
  break_pts = NULL, cor_type = "spearman", p_adj = "fdr",
  permut = 10000, ...){

  ## In sample matadata, there should be columns "Long" and "Lat"

  ## Extract coordinates
  smp <- as(sample_data(ps), "data.frame")

  ## OTU distance matrix
  odd <- phyloseq::distance(ps, method = otu_dist, type = "samples")

  ### Estimate geographic distances
  ## Euclidean distance
  if(geo_dist %in% "euclidean"){
  	gdd <- dist(smp[, c("Long", "Lat")])
  }
  ## Geodesic distance
  if(geo_dist %in% "geodesic"){
    gdd <- geodist::geodist(smp[, c("Long", "Lat")], measure = "geodesic")
    gdd <- as.dist(gdd)

    # measure = "cheap"    -- mapbox cheap ruler algo is innacurate for dists > 100 km
    # measure = "vincenty" -- errors of Vincenty distances remain constant
  }

  ## Log-transform geographic distance
  if(log_geo_dist == TRUE){ gdd <- log(gdd) }

  ## Mantel correlograms
  mant_correlog <- mantel.correlog(
    D.eco = odd,
    D.geo = gdd,
    # XY = data.xy,
    break.pts = break_pts,
    r.type = cor_type,
    cutoff = FALSE,
    nperm = permut,
    mult = p_adj, progressive = TRUE, 
    ...)

  # summary(mant_correlog)
  # plot(mant_correlog)

  ## Function to extract Mantel results
  parse_mantel <- function(x){
    # x = result of vegan::mantel.correlog
    rz <- data.frame(DistID = 1:nrow(x$mantel.res), x$mantel.res)
    rownames(rz) <- NULL
    rz$Signif <- "not_signif"
    rz$Signif[which(rz$Pr.corrected. < 0.05)] <- "signif"
    rz$Signif <-factor(rz$Signif, levels = c("not_signif", "signif"))
    return(rz)
  }

  res <- parse_mantel(mant_correlog)
  return(res)
}


# phyloseq_mantel_correlog(ps,
#   p_adj = "fdr", cor_type = "pearson", permut = 30000,
#   break_pts = log( c(1, 11, 20, 40, 80, 100, 300) )  # log breakpoints
