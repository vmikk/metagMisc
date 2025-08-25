
#' Add sister tips to a phylogenetic tree
#'
#' @description 
#' Adds one or more sister tips to an existing phylogenetic tree. 
#' Sister tips share the same immediate common ancestor with the specified existing tip(s).
#'
#' @param tree A phylogenetic tree object of class \code{phylo}
#' @param new.tip.label Character vector of labels for the new tip(s) to add
#' @param sister.tip.label Character vector of existing tip labels that will become 
#'   sisters to the new tips. If a single value is provided and multiple new tips 
#'   are specified, all new tips will be added as sisters to the same tip
#' @param branch.length Numeric value for the branch length of new tips (default: 1e-5)
#'
#' @details
#' This function supports three main scenarios: 
#' (1) adding a single tip to a single sister tip, 
#' (2) adding multiple tips to a single sister tip, and 
#' (3) adding multiple tips to multiple different sister tips. 
#' 
#' Branch lengths are automatically adjusted based on the number of tips being added to prevent very short branches 
#' that could cause numerical instability.
#'
#' @return A modified phylogenetic tree object of class \code{phylo} with the new 
#'   sister tip(s) added
#'
#' @examples
#' library(ape)
#' 
#' # Create a simple tree
#' set.seed(42)
#' tr <- rtree(n = 4)
#' 
#' # Add a single sister tip
#' tr_single <- tree_add_sister_tip(
#'   tree = tr, 
#'   new.tip.label = "NewTip", 
#'   sister.tip.label = "t4")
#' par(mfrow = c(1, 3))
#' plot(tr)
#' plot(tr_single)
#' plot(tr_single, use.edge.length = F)
#' 
#' # Add multiple tips to a single sister (creates polytomy)
#' tr_poly <- tree_add_sister_tip(
#'   tree = tr, 
#'   new.tip.label = c("Tip1", "Tip2"), 
#'   sister.tip.label = "t4")
#' par(mfrow = c(1, 3))
#' plot(tr)
#' plot(tr_poly)
#' plot(tr_poly, use.edge.length = F)
#' 
#' # Add tips to different sisters
#' tr_multi <- tree_add_sister_tip(
#'   tree = tr, 
#'   new.tip.label = c("NewTip1", "NewTip2"), 
#'   sister.tip.label = c("t4", "t3"))
#' par(mfrow = c(1, 3))
#' plot(tr)
#' plot(tr_multi)
#' plot(tr_multi, use.edge.length = F)
#'
#' @importFrom phytools bind.tip
#' @export
#' 
tree_add_sister_tip <- function(tree, new.tip.label, sister.tip.label, branch.length = 1e-5) {

    ## Validate inputs
    if (!inherits(tree, "phylo")) {
      stop("\nTree should be an object of class `phylo`\n.")
    }

    ## Convert labels to vectors
    if (!is.vector(new.tip.label)) {
      new.tip.label <- c(new.tip.label)
    }
    if (!is.vector(sister.tip.label)) {
      sister.tip.label <- c(sister.tip.label)
    }

    ## Validate inputs
    if(length(new.tip.label) < 1){
      stop("\n`new.tip.label` is not specified\n")
    }
    if(length(sister.tip.label) < 1){
      stop("\n`sister.tip.label` is not specified\n")
    }

    ## Check if vectors have equal length
    if (length(new.tip.label) > 1){
      ## If adding multiple tips to a single sister tip
      if(length(sister.tip.label) == 1){
        sister.tip.label <- rep(sister.tip.label, length(new.tip.label))
      } else if (length(sister.tip.label) != length(new.tip.label)) {
        stop("\n`new.tip.label` and `sister.tip.label` must have the same length\n")
      }
    }

    ## Calculate minimum safe branch length based on number of tips per sister
    sister_counts <- table(sister.tip.label)
    min_branch_length <- max(1e-6, branch.length)  # ensure minimum tolerance

    ## Process each pair of tips sequentially
    new_tree <- tree
    for (i in seq_along(new.tip.label)) {

        ## Find the sister tip index
        sister.index <- match(sister.tip.label[i], new_tree$tip.label)
        if (is.na(sister.index)) {
          stop(paste("\nSister tip", sister.tip.label[i], "not found in tree\n"))
        }

        ## Get the original edge length
        original_edge <- new_tree$edge.length[which(new_tree$edge[, 2] == sister.index)]
        if (is.null(original_edge)) {
          stop("\nTree must have edge lengths\n")
        }

        ## Calculate adjusted branch length for this sister tip
        n_tips <- sister_counts[sister.tip.label[i]]
        adjusted_branch_length <- min(min_branch_length, original_edge/(n_tips + 1))

        ## For subsequent tips in the same group,
        ## attach to the previous new node
        if (i > 1 && sister.tip.label[i] == sister.tip.label[i-1]) {

            ## Find the most recently added tip
            prev_tip_index <- match(new.tip.label[i-1], new_tree$tip.label)

            ## Get its parent node
            parent_node <- new_tree$edge[which(new_tree$edge[, 2] == prev_tip_index), 1]

            ## Add new tip to the parent node
            new_tree <- phytools::bind.tip(new_tree,
                                         tip.label = new.tip.label[i],
                                         where = parent_node,
                                         position = adjusted_branch_length,
                                         edge.length = adjusted_branch_length)
        } else {
            ## Standard case - create new bifurcation
            position <- original_edge - adjusted_branch_length

            new_tree <- phytools::bind.tip(new_tree,
                                         tip.label = new.tip.label[i],
                                         where = sister.index,
                                         position = position,
                                         edge.length = adjusted_branch_length)

            ## Set the sister tip branch length equal to the new tip
            sister_edge_index <- which(new_tree$edge[, 2] == match(sister.tip.label[i], new_tree$tip.label))
            new_tree$edge.length[sister_edge_index] <- adjusted_branch_length
        }
    }

    return(new_tree)
}
