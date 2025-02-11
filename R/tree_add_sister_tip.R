
## Function to add a sister tip to a tree
tree_add_sister_tip <- function(tree, new.tip.label, sister.tip.label, branch.length = 1e-5) {

    ## Validate inputs
    if (!inherits(tree, "phylo"))
        stop("\nTree should be an object of class `phylo`\n.")

    ## Find the sister tip index
    sister.index <- match(sister.tip.label, tree$tip.label)
    if (is.na(sister.index))
        stop("sister tip not found in tree")

    ## Get the original edge length to preserve total branch length
    original_edge <- tree$edge.length[which(tree$edge[, 2] == sister.index)]
    if (is.null(original_edge))
        stop("tree must have edge lengths")

    ## Calculate position to preserve total branch length
    ## We'll place the new node very close to the tip
    position <- original_edge - branch.length

    ## Add the new tip
    new_tree <- phytools::bind.tip(tree,
                        tip.label = new.tip.label,
                        where = sister.index,
                        position = position,
                        edge.length = branch.length)

    ## Set the sister tip branch length equal to the new tip
    sister_edge_index <- which(new_tree$edge[, 2] == match(sister.tip.label, new_tree$tip.label))
    new_tree$edge.length[sister_edge_index] <- branch.length

    return(new_tree)
}

## Examples
# tr <- ape::rtree(n = 4)
# plot(tr)
# plot( tree_add_sister_tip(tree = tr, new.tip.label = "TN", sister.tip.label = "t1") )
# plot( tree_add_sister_tip(tree = tr, new.tip.label = "TN", sister.tip.label = "t1"), use.edge.length = F )

