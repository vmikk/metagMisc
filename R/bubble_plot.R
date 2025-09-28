## TO DO:
##	adjust correspondance of dendrogram to table
##	add parameter to adjust ratio of dendro to table (with heights=c(x/y, z/y))

#' @title  Bubble plot
#'
#' @param x Data frame with data (columns = samples, rows = species abundance, rownames = species names)
#' @param transp Transparency level (0-1)
#' @param circ Circle scale factor
#' @param add.dendro Logical, adds dendrogram for samples
#' @param ...
#'
#' @return Invisibly returns a plot ('ggplot' class).
#' 
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous theme_bw theme element_text labs geom_segment ggplotGrob
#' @importFrom ggdendro dendro_data segment theme_dendro
#' @importFrom vegan vegdist
#' @importFrom grid unit unit.pmax
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples
#' x <- as.data.frame(matrix(runif(n = 100, min = 0, max = 100), nrow = 10))
#' rownames(x) <- paste("sp", 1:nrow(x), sep="")
#' bubble_plot(x, add.dendro=F)
#' bubble_plot(x, add.dendro=T)
#'
bubble_plot <- function(x, transp=0.9, circ=16, add.dendro = FALSE, ...){

	if(add.dendro == FALSE){
		# reshape data
		xx <- data.frame(Spec = rownames(x), x)
		setDT(xx)
		xx <- melt(data = xx, value.name = "abund", id.vars = "Spec", variable.name = "Sample")
		setDF(xx)

		# plot data
		p1 <- ggplot(xx, aes(x = Sample, y = Spec)) +
			geom_point(aes(size = abund, colour = abund), shape = 19, alpha = transp) +
			scale_size_continuous(name = "Counts ", range=c(0, circ)) +
			# scale_color_gradient(low="white", high="red") +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90)) +		# , hjust = 1
			labs(x=NULL, y=NULL)

		plot(p1)
		res <- p1
	}

	if(add.dendro == TRUE){

		# Cluster samples
		dd.row <- as.dendrogram(hclust(vegdist(t(x), method="bray")))
		row.ord <- order.dendrogram(dd.row)

		# order samples
		xx <- x[, row.ord]
		xx <- data.frame(Spec = rownames(xx), xx)
		setDT(xx)
		xx <- melt(data = xx, value.name = "abund", id.vars = "Spec", variable.name = "Sample")
		setDF(xx)

		# Extract dendrogram data and create the plots
		ddata_x <- dendro_data(dd.row)

		### Create plot components ###
		# Bubble plot
		p1 <- ggplot(xx, aes(x = Sample, y = Spec)) +
			geom_point(aes(size = abund, colour = abund), shape = 19, alpha = transp) +
			scale_size_continuous(name = "Counts ", range=c(0, circ)) +
			theme_bw() +
			theme(axis.text.x = element_text(angle = 90)) +
			labs(x=NULL, y=NULL) +
			theme(plot.margin = unit(c(0,1,1,1), "lines"))

		# Dendrogram for samples
		p2 <- ggplot(segment(ddata_x)) +
			geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
			theme_dendro() +
			theme(plot.margin = unit(c(0.2, 0, -1, 0), "lines"))

		### Draw graphic ###

		gp1 <- ggplotGrob(p1)
		gp2 <- ggplotGrob(p2)
		maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
		gp1$widths[2:5] <- as.list(maxWidth)
		gp2$widths[2:5] <- as.list(maxWidth)

		res <- grid.arrange(gp2, gp1, ncol=1, heights=c(1/5, 4/5))
	}
  invisible(res)
}
