
# ***********************************************************************************************
#   Simulation of tree distributions with varying branch length and topology
#
# ***********************************************************************************************
#   Pierre de Villemereuil & Shinichi Nakagawa (2018)
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


library(ape) 
library(MASS)
source("simul_func.R")

## ---------------- Loading tree data
load("../Trees/trees.Rdata")

## ---------------- Setting up the parameters

## Phenotypic data parameters
alpha <- 5
beta <- 2
sigma_vec <- c(2, 5, 10, 15)

## Trees distributions parameters
Ntree <- 100
brlnoise_vec <- c(0, 0.1, 0.2, 0.4, 0.7, 0.9)
nbsw_vec <- c(0, 1, 2, 5, 10, 20, 30)
R <- 100        # Number of replicates

## ---------------- Simulating the data

# For each original tree
for (t in 1:length(ltree)) {
    tree <- ltree[[t]]
    N <- length(tree[["tip.label"]])
    print(paste0("*---------------Simulations for tree ", t))
    alltips <- 1:N
    
    ## Removing the forbidden tips
    # Forbidden tips are tips that are directly connected to the root 
    # (can't swap anything with those!)
    for (i in 1:N) {
        if(dim(chain.nodes(tree, i))[1] == 1) {
            alltips <- alltips[-i]      # Removing the bad tip from the available tips
        }
    }
    
    ## Simulating tree distributions
    # For each number of swaps value
    for (nbsw in nbsw_vec) {
        print(paste0("## ---Using ", nbsw, " swaps"))
        # For each branch length noise level
        for (brlnoise in brlnoise_vec) {
            print(paste0("# Branch lenght noise:", brlnoise))
            # Set-up the swap parameters
            swap_tips <- sample(alltips, nbsw)
            swap_thres <- runif(nbsw, min = 0.1, max = 1)
            nodes <- numeric(nbsw)
            levels <- numeric(nbsw)
            
            # If we swaps, we need to check quite a few things
            if (nbsw != 0) {
                # Testing whether some nodes are equal
                OK <- FALSE
                # Continue while some thresholds are still failing
                while (!OK) {
                    # Test all nodes
                    for (k in 1:nbsw) {
                        # Finding the nodes
                        nodes[k] <- sim.swap(tree, swap_thres[k], swap_tips[k])[["node"]]
                        # While we're at it, recording the actual levels of the swaps
                        levels[k] <- sim.swap(tree, swap_thres[k], swap_tips[k])[["level"]]
                    }
                    if (length(unique(nodes)) != length(nodes)) {
                        swap_tips[duplicated(nodes)] <- sample(alltips, 
                                                               sum(duplicated(nodes)))
                        swap_thres[duplicated(nodes)] <- runif(sum(duplicated(nodes)), 
                                                               min = 0.1,
                                                               max = 1)
                    } else { OK <- TRUE }
                }
            }
            
            # Initialise the list of trees
            distree <- list()
            for (i in 1:Ntree) {
                # (Re)Initialise the tree
                phylo <- tree
                # If we need to swap
                if (nbsw != 0) {
                    # Swap nbsw times
                    for (k in 1:nbsw) {
                        # Do swap with probability 0.5
                        if (runif(1) > 0.5) {
                            phylo <- sim.swap(phylo, swap_thres[k], swap_tips[k])[["phylo"]]
                        }
                    }
                }
                # Add branch length noise
                phylo <- sim.length(phylo, brlnoise)
                # Save the result in the list of trees
                distree[[i]] <- phylo
            }
            
            # Save the list of trees (distribution)
            save(distree, 
                 file = paste0("../Distributions/trees_t", t,
                               "_nbsw", nbsw, 
                               "_brlnoise", brlnoise, ".Rdata"))
            # Save the swap info for future reference
            save(swap_thres, swap_tips, nodes, levels, 
                 file = paste0("../Distributions/swapinfo_t", t, 
                               "_nbsw", nbsw, 
                               "_brlnoise", brlnoise, ".Rdata"))
            
            ## Simulating the phenotypic data
            # For each sigma value ("residual" variance)
            for (s in sigma_vec) {
                ldata <- list()
                # Simulating R replicates
                for (r in 1:R) {
                    # Simulating "environment"
                    x <- rnorm(N, 5, 2.5)
                    # Simulating phenotypic trait according to "environment"
                    # and the original phylogenetic tree
                    y <- as.vector(mvrnorm(1, alpha + beta * x, (s^2) * cov2cor(vcv.phylo(tree))))
                    # Saving dataset for tree t and sigma s
                    ldata[[r]] <- data.frame(phen = y, env = x)
                }
                save(ldata,
                     file = paste0("../Phenotypes/data_t", t,
                                 "_s", s,
                                 "_nbsw", nbsw,
                                 "_brlnoise", brlnoise, ".Rdata"))
            }
        }
    }
}
