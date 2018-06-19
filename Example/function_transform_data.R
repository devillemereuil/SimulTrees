# ***********************************************************************************************
#   Modifications of some functions from Rphylopars
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

convert_data <- function(trait_data,tree)
{
    tr <- tree
    tr$edge.length <- tree$edge.length/mean(tree$edge.length)
    dat <- trait_data
    for(i in 2:(ncol(dat)))
    {
        dat[,i] <- (dat[,i]-mean(dat[,i],na.rm=TRUE))/sd(dat[,i],na.rm = TRUE)
    }
    return(list(trait_data=dat,tree=tr))
}
revert_data <- function(phylopars_object,original_trait_data,original_tree)
{
    phylopars_object$anc_cov <- lapply(phylopars_object$anc_cov,function(X) diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5) %*% X %*% diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5))
    phylopars_object$anc_var <- (phylopars_object$anc_var)*(matrix(1,nrow(phylopars_object$anc_var)) %*% t(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)))
    for(i in 1:ncol(phylopars_object$anc_recon)) phylopars_object$anc_recon[,i] <- phylopars_object$anc_recon[,i]*sd(original_trait_data[,i+1],na.rm=TRUE)+mean(original_trait_data[,i+1],na.rm=TRUE)
    phylopars_object$pars <- lapply(phylopars_object$pars,function(X) (diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5) %*% X %*% diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5))/mean(original_tree$edge.length))
    return(phylopars_object)
} 

## Trying using "max" instead of "mean" to scale the branch length

source("simul_func.R")

convert_data2 <- function(trait_data,tree)
{
    tr <- tree
    tree_scale <- max(sapply(1:N,function(i){sum(chain.nodes(tree,i)[,"length"])}))
    #   tr$edge.length <- tree$edge.length/mean(tree$edge.length)
    tr$edge.length <- tree$edge.length/tree_scale
    dat <- trait_data
    for(i in 2:(ncol(dat)))
    {
        dat[,i] <- (dat[,i]-mean(dat[,i],na.rm=TRUE))/sd(dat[,i],na.rm = TRUE)
    }
    return(list(trait_data=dat,tree=tr))
}
revert_data2 <- function(phylopars_object,original_trait_data,original_tree)
{
    tree_scale <- max(sapply(1:N,function(i){sum(chain.nodes(original_tree,i)[,"length"])}))
    phylopars_object$anc_cov <- lapply(phylopars_object$anc_cov,function(X) diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5) %*% X %*% diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5))
    phylopars_object$anc_var <- (phylopars_object$anc_var)*(matrix(1,nrow(phylopars_object$anc_var)) %*% t(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)))
    for(i in 1:ncol(phylopars_object$anc_recon)) phylopars_object$anc_recon[,i] <- phylopars_object$anc_recon[,i]*sd(original_trait_data[,i+1],na.rm=TRUE)+mean(original_trait_data[,i+1],na.rm=TRUE)
    #   phylopars_object$pars <- lapply(phylopars_object$pars,function(X) (diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5) %*% X %*% diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5))/mean(original_tree$edge.length))
    phylopars_object$pars <- lapply(phylopars_object$pars,function(X) (diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5) %*% X %*% diag(apply(original_trait_data[,2:ncol(original_trait_data)],2,var,na.rm=TRUE)^.5))/tree_scale)
    return(phylopars_object)
}
