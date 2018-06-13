
# ***********************************************************************************************
#   Helper functions for the rest of the script, especially for simulating the trees
#
# ***********************************************************************************************
#                       Pierre de Villemereuil & Shinichi Nakagawa (2018)
# ***********************************************************************************************
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -------------------------------- Functions for simulation

## Modifying the branch length by a certain factor (from 0 to 1)
sim.length <- function(phylo, factor) {
    # Check that factor is OK
    if (factor < 0 | factor > 1) { stop("Factor should be between 0 and 1") }
    
    # Add noise to branch length
    noise <- factor * phylo[["edge.length"]]
    phylo[["edge.length"]] <- phylo[["edge.length"]] + runif(length(noise), -noise, noise)
    return(phylo)
}

## Finding the chain of nodes from tip to root
chain.nodes <- function(phylo, tip) {
    # Check that tip is OK
    if (tip > length(phylo[["tip.label"]])) { stop("Tip is too large!") }
    
    # Setting the parameters up
    continue  <- TRUE     # Boolean for the while lopp
    node      <- tip      # Tip is the first node
    length    <- NULL     # will contain the length of all edges
    all       <- NULL     # will contain all nodes
    
    # Get the nodes up to the root
    while (continue) {
        # Stocking nodes
        all <- c(all, node)
        # Length of the nodes
        length <- c(length, phylo[["edge.length"]][which.edge(phylo, node)])
        # Find the next node (which.edge looks at the second column of edge for "node")
        node <- phylo[["edge"]][which.edge(phylo, node), 1]
        # If the next node is not found by which.edge, we found the root, so we stop
        if (is.na(which.edge(phylo, node))) { continue <- FALSE }
    }
    
    # Return a data.frame with nodes, length of the edges 
    # and cumulative length from tip to root (standardised to be of one)
    # Caution: the cumulative length will always go up to 1, even if the tree is not ultrametric!
    data.frame(node = all, length = length, cumlength = cumsum(length)/sum(length, na.rm = T))
}

## Finding the chain of nodes from a node to tip
chain.nodes.reverse <- function(phylo, nodes) {
    # Setting the parameters up
    continue <- TRUE      # Boolean for the while loop
    Nspec <- length(phylo[["tip.label"]])
    all <- NULL           # will contain all tips
    
    # Finding the chain of nodes
    while (continue) {
        # Find the next node
        nodes <- phylo[["edge"]][phylo[["edge"]][, 1] %in% nodes, 2]
        if (sum(nodes < Nspec) > 0) {
            all <- c(all, nodes[nodes < Nspec])
            nodes <- nodes[nodes > Nspec]
        }
        if (sum(nodes > Nspec) == 0) {
            continue <- FALSE
        }
    }
    # Return tips
    return(all)
}

## Function swapping some nodes in the tree
sim.swap <- function(phylo, threshold, tip = NULL) {
    # Check that threshold is OK
    if (threshold < 0 | threshold > 1) { stop("Threshold should be between 0 and 1") }
    
    # If no tip provided, randomly select a tip
    if (is.null(tip)) { tip <- sample(1:length(phylo[["tip.label"]]), 1) }
    
    # Check that tip is OK
    if (tip > length(phylo[["tip.label"]])) { stop("Tip is too large!") }
    
    # Retrace the path from tip to root (see chain.nodes above)
    chain <- chain.nodes(phylo, tip)
    # Selecting the most upper node right below the treshold
    node <- tail(chain[chain[["cumlength"]] < threshold, "node"], 1)
    # Finding the node just below "node" in the path
    previousnode <- tail(chain[chain[["cumlength"]] < threshold, "node"], 2)[1]
    
    # If threshold is way too low, then node is empty
    if (length(node) == 0) {
        node <- tip		# Dirty fix, simply allow to go into the second particular case
    }
    
    # We have a problem when node is equal to tip
    if (node == tip) {
        # Node is set to be the "next node in line"
        node <- phylo[["edge"]][which.edge(phylo, node), 1]
        # And previousnode is now tip (normally, it should already be tip, but well...)
        previousnode <- tip
    }
    
    # Finding the node just above "node" in the path
    nextnode <- phylo[["edge"]][which.edge(phylo, node), 1]
    # "nextnode" is branching in at least two edges -> one goes to "node", 
    # the other(s) are going to sister group(s)
    vecgroup <- phylo[["edge"]][phylo[["edge"]][, 1] == nextnode, 2]
    
    # removing "node" from there
    vecgroup <- vecgroup[vecgroup != node]
    
    # Selecting a sister group
    sistergroup <- vecgroup[1]
    if (sistergroup <= length(phylo[["tip.label"]])) {
        sistergroup <- phylo[["edge"]][which.edge(phylo, sistergroup), 1]
    }
    
    # Swapping the edges: 
    # breaking the edge between "previousnode" and "node" 
    # and creating an edge between "previousnode" and "sistergroup"
    phylo[["edge"]][phylo[["edge"]][, 2] == previousnode, 1] <- sistergroup
    
    # We created a non-branching node, which ape doesn't like, fixing that:
    phylo <- collapse.singles(phylo)
    # Recording the actual level of the swap
    level <- chain[chain[["node"]] == previousnode, "cumlength"]
    
    # Return the transformed phylogeny
    return(list(phylo = phylo, level = level, node = node))
} 

# find.swap <- function(phylo, threshold, tip) {
#   # Check that threshold is OK
#   if (threshold < 0 | threshold > 1) { stop("Threshold should be between 0 and 1") }
#
#   # If not tip provided, randomly select a tip
#   if (is.null(tip)) { tip <- sample(1:length(phylo[["tip.label"]]), 1) }
#
#   # Check that tip is OK
#   if (tip > length(phylo[["tip.label"]])) { stop("Tip is too large!") }
#
#   # Retrace the path from tip to root (see chain.nodes above)
#   chain <- chain.nodes(phylo, tip)
#
#   # Selecting the most upper node right below the treshold
#   node <- tail(chain[chain[["cumlength"]] < threshold, "node"], 1)
#   return(node)
# }
