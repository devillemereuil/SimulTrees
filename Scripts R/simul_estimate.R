
# ***********************************************************************************************
#   Simulation of datasets to study the distribution of estimates according to the method
#
# ***********************************************************************************************
#                   Pierre de Villemereuil & Shinichi Nakagawa (2018)
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
#    along with this program.  If not, see <http://www.gnu.org/licenses/ > .

library(ape)
library(mice)
library(nlme)
library(AICcmodavg)

##---------------Setting up parameters
# going through - 
sigmaVec        <- c(2, 5, 10, 15)
which_treeVec   <- 1:13
nbswVec         <- c(0, 1, 2, 5, 10, 20, 30)
brlnoiseVec     <- c(0, 0.1, 0.2, 0.4, 0.7, 0.9)

for(h in 1:4){
    sigma <- sigmaVec[h]
    print(paste0("sigma", sigma))
    
    for(i in 1:13){
        which_tree <- which_treeVec[i]
        print(paste0("tree", which_tree))
        
        for(j in 1:7){
            nbsw <- nbswVec[j]
            print(paste0("swap", nbsw))
            
            for(k in 1:6){
                brlnoise <- brlnoiseVec[k]
                print(paste0("noise", brlnoise))
                
                #----------------------
                #Fixed parameters
                #Number of trees
                M <- 100
                #Number of replicates
                R <- 100
                #Regression parameters
                alpha <- 5
                beta  <- 2 
                
                ##---------------Loading the data
                #Loading the true tree
                load("Trees/trees.Rdata")
                true_tree <- ltree[[which_tree]]
                true_mat <- cov2cor(vcv.phylo(true_tree))
                
                #Loading the distribution of trees
                load(paste0("Distributions/trees_t", which_tree,
                            "_nbsw", nbsw, 
                            "_brlnoise", brlnoise, ".Rdata"))
                distmat <- lapply(distree, function(tree) { cov2cor(vcv.phylo(tree)) })
                #Loading the phenotypic data
                load(paste0("Phenotypes/data_t", which_tree,
                            "_s", sigma,
                            "_nbsw", nbsw,
                            "_brlnoise", brlnoise, ".Rdata"))
                
                #Number of species
                N <- length(true_tree[["tip.label"]])
                
                #Constructing the consensus tree
                constree_strict <- compute.brlen(consensus(distree))
                consmat_strict  <- cov2cor(vcv.phylo(constree_strict))
                constree_maj    <- compute.brlen(consensus(distree, p = 0.5))
                consmat_maj     <- vcv.phylo(constree_maj, corr = TRUE)
                
                ##----------------Performing the analysis
                
                for (r in 1:R) {
                    
                    data <- ldata[[r]]
                    
                    ##Fitting simple GLSs
                    #True model
                    true_model <- gls(phen ~ env,
                                      correlation = 
                                          corSymm(true_mat[lower.tri(true_mat)],
                                                  fixed = TRUE),
                                      data = data)
                    
                    #Strict consensus tree
                    strictcons_model <- gls(phen ~ env, 
                                            correlation =
                                                corSymm(consmat_strict[lower.tri(consmat_strict)],
                                                        fixed = TRUE),
                                            data = data)
                    
                    #Majority rule consensus tree
                    majcons_model <- gls(phen ~ env,
                                         correlation = 
                                             corSymm(consmat_maj[lower.tri(consmat_maj)],
                                                     fixed = TRUE),
                                         data = data)
                    
                    ##Fitting multiple GLS
                    #Fitting the models
                    mods <- list()
                    for (m in 1:M) {
                        mods[[m]] <- gls(phen ~ env,
                                         correlation = 
                                             corSymm(distmat[[m]][lower.tri(distmat[[m]])],
                                                     fixed = TRUE),
                                         data = data)
                    }
                    
                    #Combining estimates following Rubin's rule
                    models <- as.mira(mods)
                    pool_small <- summary(pool(models, method = "smallsample"))
                    pool_orig <- summary(pool(models, method = "original"))
                    
                    ##AIC averaging
                    AIC <- sapply(mods, function(model){
                        -2 * summary(model)[["logLik"]] + 2 * 2 * N / (N - 2 - 1)
                    })
                    deltaAIC    <- AIC - min(AIC)
                    weightAIC   <- exp(-deltaAIC / 2) / sum(exp(-deltaAIC / 2))
                    avg_alpha   <- modavg(parm = "(Intercept)", cand.set = mods)
                    avg_beta    <- modavg(parm = "env", cand.set = mods)
                    
                    #Coefficients
                    coeffs_true <- 
                        c(true_model[["coefficients"]],
                          sigma = true_model[["sigma"]])
                    coeffs_strictcons <- 
                        c(strictcons_model[["coefficients"]],
                          sigma = strictcons_model[["sigma"]])
                    coeffs_majcons <- 
                        c(majcons_model[["coefficients"]],
                          sigma = majcons_model[["sigma"]])
                    coeffs_multi_small <- 
                        c(pool_small[ , "est"], 
                          sigma = mean(sapply(mods, function(model) {
                              summary(model)[["sigma"]]
                          })))
                    coeffs_multi_orig <- 
                        c(pool_orig[ , "est"],
                          sigma = mean(sapply(mods, function(model) {
                              summary(model)[["sigma"]]
                          })))
                    coeffs_AIC <- 
                        c(avg_alpha[["Mod.avg.beta"]],
                          avg_beta[["Mod.avg.beta"]],
                          sum(sapply(mods, 
                                     function(model) { summary(model)[["sigma"]] }) * weightAIC))
                    
                    
                    #Confidence Interval
                    upper_true <- intervals(true_model)[["coef"]][ , 'upper']
                    lower_true <- intervals(true_model)[["coef"]][ , 'lower']
                    upper_strictcons <- intervals(strictcons_model)[["coef"]][ , 'upper']
                    lower_strictcons <- intervals(strictcons_model)[["coef"]][ , 'lower']
                    upper_majcons <- intervals(majcons_model)[["coef"]][ , 'upper']
                    lower_majcons <- intervals(majcons_model)[["coef"]][ , 'lower']
                    upper_multi_small <- pool_small[ , "hi 95"]
                    lower_multi_small <- pool_small[ , "lo 95"]
                    upper_multi_orig <- pool_orig[ , "hi 95"]
                    lower_multi_orig <- pool_orig[ , "lo 95"]
                    upper_AIC <- c(avg_alpha[["Upper.CL"]], avg_beta[["Upper.CL"]])
                    lower_AIC <- c(avg_alpha[["Lower.CL"]], avg_beta[["Lower.CL"]])
                    
                    #Tests
                    bool_true <- 
                        c((upper_true[1] < alpha) | (lower_true[1] > alpha),
                          (upper_true[2] < beta)  | (lower_true[2] > beta))
                    bool_strictcons <- 
                        c((upper_strictcons[1] < alpha) | (lower_strictcons[1] > alpha),
                          (upper_strictcons[2]   < beta)  | (lower_strictcons[2] > beta))
                    bool_majcons <- 
                        c((upper_majcons[1] < alpha) | (lower_majcons[1] > alpha),
                          (upper_majcons[2] < beta)  | (lower_majcons[2] > beta))
                    bool_multi_small <- 
                        c((upper_multi_small[1] < alpha) | (lower_multi_small[1] > alpha),
                          (upper_multi_small[2]   < beta)  | (lower_multi_small[2] > beta))
                    bool_multi_orig <- 
                        c((upper_multi_orig[1] < alpha) | (lower_multi_orig[1] > alpha),
                          (upper_multi_orig[2]   < beta)  | (lower_multi_orig[2] > beta))
                    bool_AIC <- 
                        c((upper_AIC[1] < alpha) | (lower_AIC[1] > alpha),
                          (upper_AIC[2] < beta)  | (lower_AIC[2] > beta))
                    
                    #Saving estimates
                    write(coeffs_true, 
                          file = paste0("./Output/trueGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    write(coeffs_strictcons, 
                          file = paste0("./Output/strictconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    write(coeffs_majcons,
                          file = paste0("./Output/majconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    write(coeffs_multi_small, 
                          file = paste0("./Output/multiGLS_small_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    write(coeffs_multi_orig, 
                          file = paste0("./Output/multiGLS_orig_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    write(coeffs_AIC,
                          file = paste0("./Output/avgGLS_AIC_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".estimate"),
                          append = TRUE)
                    
                    #Saving CI
                    write(c(upper_true, lower_true),
                          file = paste0("./Output/trueGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    write(c(upper_strictcons, lower_strictcons),
                          file = paste0("./Output/strictconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    write(c(upper_majcons, lower_majcons),
                          file = paste0("./Output/majconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    write(c(upper_multi_small, lower_multi_small),
                          file = paste0("./Output/multiGLS_small_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    write(c(upper_multi_orig, lower_multi_orig),
                          file = paste0("./Output/multiGLS_orig_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    write(c(upper_AIC, lower_AIC),
                          file = paste0("./Output/avgGLS_AIC_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".interv"),
                          append = TRUE,
                          ncolumns = 6)
                    
                    #Saving tests
                    write(as.numeric(bool_true),
                          file = paste0("./Output/trueGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                    write(as.numeric(bool_strictcons),
                          file = paste0("./Output/strictconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                    write(as.numeric(bool_majcons),
                          file = paste0("./Output/majconsGLS_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                    write(as.numeric(bool_multi_small),
                          file = paste0("./Output/multiGLS_small_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                    write(as.numeric(bool_multi_orig),
                          file = paste0("./Output/multiGLS_orig_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                    write(as.numeric(bool_AIC),
                          file = paste0("./Output/avgGLS_AIC_t", which_tree,
                                        "_nbsw", nbsw,
                                        "_brlnoise", brlnoise,
                                        "_s", sigma, ".test"),
                          append = TRUE)
                }
            }
        }
    }
}
