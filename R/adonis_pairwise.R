
# x = original data
# dd = global distance matrix
# group.var = name of the groupping variable
# p.adj = Adjust P-values for Multiple Comparisons?
# adj.meth = "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"

adonis_pairwise <- function(x, dd, group.var = "Fact", permut = 999,
	p.adj=T, adj.meth="fdr", ...){
	# expand.grid(levels(x[,1]), levels(x[,1]))

	VV <- which(colnames(x) %in% group.var)
	cc <- combn(levels(x[,VV]), 2)					# combinations to test

	dd <- as.matrix(dd)

	# subset distance matrix
	prep.dist <- function(x, dd, sub1, sub2){
		subs <- which(x[,VV] == sub1 | x[,VV] == sub2)
		gr <- x[,VV][subs]
		res <- list()
			res$dd <- as.dist(dd[subs, subs])
			res$gr <- gr
		return(res)
	}

	# prepare data for pair-wise comparisons
	dd.subs <- list()
	dd.groups <- list()
	for(i in 1:ncol(cc)){
		tmp <- prep.dist(x, dd, sub1=cc[1,i], sub2=cc[2,i])
		dd.subs[[i]] <- tmp$dd
		dd.groups[[i]] <- tmp$gr
		rm(tmp)
	}

	# pair-wise adonis
	adon <- list()
	for(i in 1:length(dd.subs)){
		adon[[i]] <- adonis(dd.subs[[i]] ~ dd.groups[[i]], permutations=permut, ...)
	}

	# prepare names
	foo <- vector()
	for(i in 1:ncol(cc)){
		foo[i] <- paste(cc[,i], collapse = ".")
	}
	names(dd.subs) <- names(dd.groups) <- names(adon) <- foo

	# extract results
	adon.extract.p <- function(adon){ adon$aov.tab$Pr[1] }
	adon.extract.F <- function(adon){ adon$aov.tab$F.Model[1] }
	adon.extract.df <- function(adon){ paste(adon$aov.tab$Df[1:2], collapse=";") }
	ad.t <- data.frame(Comparison = foo,
					F = unlist(lapply(adon, FUN = adon.extract.F)),
					df = unlist(lapply(adon, FUN = adon.extract.df)),
					p = unlist(lapply(adon, FUN = adon.extract.p)))
	rownames(ad.t) <- NULL

	if(p.adj == TRUE){
		ad.t$p.adj <- p.adjust(ad.t$p, method = adj.meth)
	}

	res <- list()
		res$Adonis.tab <- ad.t
		res$Adonis <- adon
		res$Dist.subsets <- dd.subs
		res$Groups <- dd.groups
	return(res)
}

## Example
library(vegan)
data(dune)
data(dune.env)

# Compare all Management levels
adonis(dune ~ Management, data = dune.env)

# Pairwise comparisons between Management levels
ad <- adonis_pairwise(x = dune.env, dd = vegdist(dune), group.var = "Management")
ad$Adonis.tab

