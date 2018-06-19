
# ***********************************************************************************************
#   Example file illustrating how to use Rubin's rule to perform phylogenetic comparative
#                                       analysis in R
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

## ----------- Loading the libraries needed
library(ape)
library(nlme)
library(mice)
library(Rphylopars)

# Using older version of pool when mice 3.0 is used
source("old.pool.mice.R")


## ********************************************************************************************
##                                      Preliminary steps
## ********************************************************************************************

## ----------- Loading the example dataset

# Loading the tree distribution
load("trees.Rdata")
# This is tree #11, nb swaps 10 and br. lenth noise of 0.2
# The object distree contains a list of 100 trees
head(distree)

# Loading the data
load("data.Rdata")
# Data contains the phenotypic values and the environmental covariate
head(data)

## ----------- Preparing the analysis

# Computing the VCV phylogenetic matrices
distmat <- lapply(distree, function(tree) { cov2cor(vcv.phylo(tree)) })
# Note that we use correlation matrices to obtain a variance relating to the variance of the data

## ********************************************************************************************
##                  Accounting only for phylogenetic uncertainty
## ********************************************************************************************

## ----------- Performing the analysis

# Fitting the models
mods <- list()
for (m in 1:length(distree)) {
    mods[[m]] <- gls(phen ~ env,
                     # Using the VCV phylogenetic matrix
                     correlation = corSymm(distmat[[m]][lower.tri(distmat[[m]])], fixed = TRUE),
                     data = data)
}

# Combining estimates following Rubin's rule
models <- as.mira(mods)
# as.mira takes the list of models and create an object to be used by the mice package
pool   <- summary(pool(models, method = "smallsample"))
# pool summarise the models using Rubin's rule corrected for small samples

# We now can look at the estimates
pool
# For each parameter, we have the estimate values (est) and their standard errors (se),
# we also have a significance test (t, df and  Pr(>|t|)) and 95% CI (lo 95, hi 95)
# based on the degrees of freedom computed using Rubin's rule.
# We then have information about the multiple imputation such as number of missing values (nmis),
# fraction of missing information (fmi) and lambda
# (for more information, see the article this example file is linked to)

# We can compute the efficiency of the multiple imputation as
eff <- 1 / (1 + (max(pool[ , "fmi"]) / length(distree)))
# Should be pretty close to 1

## ********************************************************************************************
##                          Accounting for missing data as well
## ********************************************************************************************
    
# Simulate missing data
p <- 0.1
N <- nrow(data)
Nmiss <- round(p * N)
data[["phen"]][sample(1:N, Nmiss)] <- NA

# Use functions to convert data
source("function_transform_data.R")

# Looping over the trees
mods <- list()
for (l in 1:length(distree)) {
    ##TODO Finish this and find a way to make Rphylopars comply!
    dat  <- convert_data2(trait_data = data,
                          tree = distree[[l]])
    pvcv <- revert_data2(phylopars_object = 
                             phylopars(trait_data = dat[["trait_data"]],
                                       tree = dat[["tree"]]),
                         original_trait_data = data,
                         original_tree = distree[[l]])
}
