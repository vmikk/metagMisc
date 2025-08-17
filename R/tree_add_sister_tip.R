
## Function to add a sister tip to a tree
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
                                         position = branch.length,
                                         edge.length = branch.length)
        } else {
            ## Standard case - create new bifurcation
            position <- original_edge - branch.length

            new_tree <- phytools::bind.tip(new_tree,
                                         tip.label = new.tip.label[i],
                                         where = sister.index,
                                         position = position,
                                         edge.length = branch.length)

            ## Set the sister tip branch length equal to the new tip
            sister_edge_index <- which(new_tree$edge[, 2] == match(sister.tip.label[i], new_tree$tip.label))
            new_tree$edge.length[sister_edge_index] <- branch.length
        }
    }

    return(new_tree)
}

set.seed(42)
tr <- ape::rtree(n = 4)

## Add a single tip
trs <- tree_add_sister_tip(tree = tr, new.tip.label = "TN", sister.tip.label = "t4")
par(mfrow = c(1, 3))
plot(tr)
plot(trs)
plot(trs, use.edge.length = F)

## Add multiple tips to a single sister tip
trs <- tree_add_sister_tip(tree = tr, new.tip.label = c("TN1", "TN2"), sister.tip.label = "t4")
par(mfrow = c(1, 3))
plot(tr)
plot(trs)
plot(trs, use.edge.length = F)

## Independently add several tips to several sister tips
trs <- tree_add_sister_tip(tree = tr, new.tip.label = c("TN4", "TN3"), sister.tip.label = c("t4", "t3"))
par(mfrow = c(1, 3))
plot(tr)
plot(trs)
plot(trs, use.edge.length = F)

