
# ***********************************************************************************************
#   Simulation of datasets to study the efficiency of multiple imputation
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

library(ape)
library(mice)
library(nlme)
library(AICcmodavg)

## ---------------Setting up parameters
# going through - 
sigmaVec        <- c(2, 5, 10, 15)
which_treeVec   <- 2:13
nbswVec         <- c(0, 1, 2, 5, 10, 20, 30)
brlnoiseVec     <- c(0, 0.1, 0.2, 0.4, 0.7, 0.9)
nbtreeVec       <- c(10, 20, 50, 100)

# Fixed parameters
# Number of replicates
R <- 100
# Regression parameters
alpha <- 5
beta  <- 2 

# Initialising the variable stocking the results
Ntot <-length(sigmaVec) *
       length(which_treeVec) *
       length(nbswVec) *
       length(brlnoiseVec) *
       length(nbtreeVec) *
       100
eff <- data.frame(sigma     = rep(NA, Ntot),
                  tree      = rep(NA, Ntot),
                  nbsw      = rep(NA, Ntot),
                  brlnoise  = rep(NA, Ntot),
                  nbtree    = rep(NA, Ntot),
                  fmi       = rep(NA, Ntot),
                  lambda    = rep(NA, Ntot))
index <- 1  # index will increase over all the loops 

for(t in 1:4) {
    nbtree <- nbtreeVec[t]
    print(paste0("NbTree", nbtree))
    
    for(h in 1:4) {
        sigma <- sigmaVec[h]
        print(paste0("sigma", sigma))
        
        for(i in 1:12) {
            which_tree <- which_treeVec[i]
            print(paste0("tree", which_tree))
            
            for(j in 1:7) {
                nbsw <- nbswVec[j]
                print(paste0("swap", nbsw))
                
                for(k in 1:6) {
                    brlnoise <- brlnoiseVec[k]
                    print(paste0("noise", brlnoise))
                    
                    
                    ## ---------------Loading the data
                    # Loading the true tree
                    load("../Trees/trees.Rdata")
                    true_tree <- ltree[[which_tree]]
                    true_mat  <- cov2cor(vcv.phylo(true_tree))
                    # Loading the distribution of trees
                    load(paste0("../Distributions/trees_t", which_tree,
                                "_nbsw", nbsw,
                                "_brlnoise", brlnoise, ".Rdata"))
                    distree <- distree[sample(1:100, nbtree)]
                    distmat <- lapply(distree, function(tree) { cov2cor(vcv.phylo(tree)) })
                    # Loading the phenotypic data
                    load(paste0("../Phenotypes/data_t", which_tree,
                                "_s", sigma,
                                "_nbsw", nbsw,
                                "_brlnoise", brlnoise, ".Rdata"))
                    
                    # Number of species
                    N <- length(true_tree[["tip.label"]])
                    
                    # Number of trees
                    M <- nbtree
                    
                    ## ----------------Performing the analysis
                    
                    for (r in 1:R) {
                        
                        data <- ldata[[r]]
                        
                        ## Fitting multiple GLS
                        # Fitting the models
                        mods <- list()
                        for (m in 1:M) {
                            mods[[m]] <- gls(phen~env, 
                                             correlation = 
                                                 corSymm(distmat[[m]][lower.tri(distmat[[m]])],
                                                         fixed = TRUE), 
                                             data = data)
                        }
                        
                        # Combining estimates following Rubin's rule
                        pool_small <- summary(pool(as.mira(mods), method = "smallsample"))
                        
                        # Saving efficiency
                        eff[index, ] <- data.frame(sigma    = sigma,
                                                   tree     = which_tree,
                                                   nbsw     = nbsw,
                                                   brlnoise = brlnoise,
                                                   nbtree   = nbtree,
                                                   fmi      = pool_small["env", "fmi"],
                                                   lambda   = pool_small["env", "lambda"])
                        index <- index + 1
                        
                        }
                }
            }
        }
    }
}

save(eff, file = "../efficiency.Rdata")
