
# ***********************************************************************************************
#   Simulation of datasets with varying missing data proportion and processes
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
library(mvtnorm)
# install.packages("Rphylopars")
library(Rphylopars) 
# install.packages("splancs") # you need this for the PVR package
# install.packages("~/Documents/PVR_0.2.1.tar.gz", repos = NULL, type = "source") # you need this from the web
library(PVR)
library(parallel)
# options(mc.cores = 8)   #set the number of parallel processes
source("simul_func.R")
source("function_transform_data.R")

##---------------Setting up parameters
# going through
sigma <- 5
which_treeVec <-2:13
nbsw <- 2
brlnoise <- 0.2
nbtree <- 100
propMissVec <- c(0.1, 0.3, 0.5)

vecmeth <- c("MCAR", "MARvar", "MARphylo")

#Fixed parameters
# Number of multiple imputation
M = 100
#Number of replicates
R <- 100
#Regression parameters
alpha <- 5
beta  <- 2 

#Initialising the variable stocking the results
Ntot <- length(which_treeVec) * length(propMissVec) * length(vecmeth) * R * 2
missing <- data.frame(tree      = rep(NA, Ntot),
                      mechanism = rep(NA, Ntot),
                      propMiss  = rep(NA, Ntot),
                      method    = rep(NA, Ntot),
                      issue     = rep(NA, Ntot),
                      intercept = rep(NA, Ntot),
                      min.intercept = rep(NA, Ntot),
                      max.intercept = rep(NA, Ntot),
                      fmi.intercept = rep(NA, Ntot),
                      slope     = rep(NA, Ntot),
                      min.slope = rep(NA, Ntot),
                      max.slope = rep(NA, Ntot),
                      fmi.slope = rep(NA, Ntot))
index <- 1   #index will increase over all the loops 

##----------------Starting the loops over the trees

for (i in 1:12) {
    
    which_tree <- which_treeVec[i]
    print(paste0("tree", which_tree))

    ##---------------Loading the data
    #Loading the true tree
    load("../Trees/trees.Rdata")
    true_tree <- ltree[[which_tree]]
    true_mat  <- cov2cor(vcv.phylo(true_tree))
    #Loading the distribution of trees
    load(paste0("../Distributions/trees_t", which_tree,
                "_nbsw", nbsw,
                "_brlnoise", brlnoise, ".Rdata"))
    distmat <- lapply(distree, function(tree) { cov2cor(vcv.phylo(tree)) })
    #Loading the phenotypic data
    load(paste0("../Phenotypes/data_t", which_tree,
                "_s", sigma,
                "_nbsw", nbsw,
                "_brlnoise", brlnoise, ".Rdata"))
    
    #Number of species
    N <- length(true_tree[["tip.label"]])
    
    #Number of trees
    L <- nbtree 
    
    ##----------------Starting the loops over datasets
    
    for (r in 1:R) {
        
        for (p in propMissVec) {
            
            for (mech in vecmeth) {
                
                #Loading the replicate dataset
                data <- ldata[[r]]
                data[ , "species"] <- true_tree[["tip.label"]]
                data <- data[c("species", "phen", "env")]
                data1 <- data                                    # need this later
                
                ##----------------Creating a missing value dataset
                
                if (mech == "MCAR") {
                    Nmiss <- round(p * N)
                    data[["phen"]][sample(1:N, Nmiss)] <- NA
                }
                
                if (mech == "MARvar") {
                    Nmiss <- round(p * N)
                    weight <- rank(data[["env"]])/nrow(data)
                    data[["phen"]][sample(x = 1:N, size = Nmiss, prob = weight)] <- NA
                }
                
                if (mech == "MARphylo") {
                    lat <- rmvnorm(1, rep(0, N), true_mat)[1, ]
                    thres <- quantile(lat, probs = p)
                    data[["phen"]][lat < thres] <- NA
                }
                
                
                ##----------------Performing the analysis
                
                ## **************** Using phylopars
                #Fitting the models
                
                print("Using phylopars")
                
                mods <- do.call("c", mclapply(1:L, function(l) {
                    
                    #looping over a few number of trials
                    good  <- FALSE
                    trial <- 1
                    while (!good & trial <= 10) {
                        # getting missing data and their variances
                        dat  <- convert_data2(trait_data = data,
                                              tree = distree[[l]])
                        pvcv <- revert_data2(phylopars_object = 
                                                 phylopars(trait_data = dat[["trait_data"]],
                                                           tree = dat[["tree"]]),
                                             original_trait_data = data,
                                             original_tree = distree[[l]])

                        #checking if phylopars ran correctly
                        if (sum(pvcv[["anc_var"]][1:dim(data)[1], 1] < 0) == 0) {
                            good <- TRUE
                        } else {trial <- trial +1}
                    }
                    
                    if (good) {
                        # get 100 sets of recovered phen using phylopars output
                        rep_phen <- 
                            replicate(M, rnorm(dim(data)[1], 
                                               mean = pvcv[["anc_recon"]][1:dim(data)[1], 1],
                                               sd = sqrt(pvcv[["anc_var"]][1:dim(data)[1], 1])))
                        
                        # using these phen values to make 100 datasets
                        rep_data <- lapply(1:M, function(m) {
                            data1[["phen"]] <- rep_phen[, m]
                            return(data1)
                        })
                        
                        #fitting gls models - should creat 10, 000 results
                        tmp <- lapply(rep_data, function(df) {
                            gls(phen ~ env,
                                correlation = corSymm(distmat[[l]][lower.tri(distmat[[l]])],
                                                      fixed = TRUE),
                                data = df)
                        })
                        if (sum(sapply(tmp, function(mod) {class(mod)}) == "gls") != length(tmp)) {
                            out_context <<- list(tree = which_tree,
                                                 mechanism = mech,
                                                 rep = r,
                                                 propMiss = p,
                                                 l = l)
                            out_pvcv <<- pvcv
                            out_rphen <<- rep_phen
                            out_rdata <<- rep_data
                            models <<- tmp
                            save(out_context, 
                                 out_pvcv,
                                 out_rphen,
                                 out_rdata,
                                 models,
                                 file = "dump.Rdata")
                            stop(paste("Error using l =", l))
                        }
                    } else {
                        tmp <- NULL
                    }
                    return(tmp)
                }))
                
                
                if (is.null(mods)) {
                    missing[index, ] <- data.frame(
                        tree        = which_tree,
                        mechanism   = mech,
                        propMiss    = p,
                        method      = "PhyloPars",
                        issue       = TRUE,
                        intercept   = NA,
                        min.intercept = NA,
                        max.intercept = NA,
                        fmi.intercept = NA,
                        slope       = NA,
                        min.slope   = NA,
                        max.slope   = NA,
                        fmi.slope   = NA)
                } else {
                    models <- as.mira(mods)
                    results1 <- summary(pool(models, method = "original"))
                    
                    issue <- ifelse(length(mods) == 10000, FALSE, TRUE)
                    
                    #Saving results
                    missing[index, ] <- data.frame(
                        tree        = which_tree,
                        mechanism   = mech,
                        propMiss    = p,
                        method      = "PhyloPars",
                        issue       = issue,
                        intercept   = results1["(Intercept)", "est"],
                        min.intercept = results1["(Intercept)", "lo 95"],
                        max.intercept = results1["(Intercept)", "hi 95"],
                        fmi.intercept = results1["(Intercept)", "fmi"],
                        slope       = results1["env", "est"],
                        min.slope   = results1["env", "lo 95"],
                        max.slope   = results1["env", "hi 95"],
                        fmi.slope   = results1["env", "fmi"]
                    )
                }
                print(paste0(index, "/", Ntot))
                print(missing[index, ])
                index <- index + 1
                
                gc()
                
                
                ## **************** Using MICE
                #Fitting the models
                print("Using MICE")
                
                # this loop takes long time as we are MIing 100 times for tree
                mods <- do.call("c", mclapply(1:L, function(l) {
                    # first get eigen vectors from a tree
                    eigen <- PVRdecomp(distree[[l]], scale = TRUE)
                    
                    eigen_vec <- as.data.frame(eigen@Eigen[["vectors"]][, 1:10])
                    eigen_vec[["species"]] <- eigen@phylo[["tip.label"]]
                    eigen_data <- merge(data, eigen_vec, by.x = "species", by.y = "species")
                    
                    # multiple imputation creating 100 data sets (per one tree)
                    imp <- mice(eigen_data[, -1], m = M, printFlag = F)
                    # analysing these 100 datasets
                    analysis <- with(imp, 
                                     gls(phen ~ env,
                                         correlation = corSymm(distmat[[l]][lower.tri(distmat[[l]])],
                                                               fixed = TRUE)))
                    # we should again get 10, 000 model outputs
                    return(analysis[[4]]) # getting model outputs
                }))
                
                
                models <- as.mira(mods)
                
                results2 <- summary(pool(models, method = "original"))
                
                #Saving results
                missing[index, ] <- data.frame(
                    tree        = which_tree,
                    mechanism   = mech,
                    propMiss    = p,
                    method      = "Eigen",
                    issue       = NA,
                    intercept   = results2["(Intercept)", "est"],
                    min.intercept = results2["(Intercept)", "lo 95"],
                    max.intercept = results2["(Intercept)", "hi 95"],
                    fmi.intercept = results2["(Intercept)", "fmi"],
                    slope       = results2["env", "est"],
                    min.slope   = results2["env", "lo 95"],
                    max.slope   = results2["env", "hi 95"],
                    fmi.slope   = results2["env", "fmi"])
                print(paste0(index, "/", Ntot))
                print(missing[index, ])
                index <- index + 1
                
                gc()
            }
        }
    }
    
    save(missing, file = "../missing.Rdata")
    
}

