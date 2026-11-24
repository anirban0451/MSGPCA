## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE,comment='#>')

## ----fit----------------------------------------------------------------------
library(MSGPCA)
sim <- simulation_gg_blocks(nside=6,Jsp=12,Jmix=24,Jns=12,seed=101)
Y <- lapply(sim$Y_list,log1p)
fit <- MSGPCA(Y,sim$locs_list,n_shared=5,n_local=2,
              lambda_shared=.3,lambda_local=.1,seed=101,verbose=FALSE)
summary(fit)
plot_convergence(fit)

## ----clustering---------------------------------------------------------------
gmm <- cluster_MSGPCA(fit,5,method='gmm',seed=101)
km <- cluster_MSGPCA(fit,5,method='kmeans',seed=101)
plot_clusters(gmm)
head(gmm$raw[[1]])
head(gmm$refined[[1]])

## ----downstream---------------------------------------------------------------
shared <- reconstruct_MSGPCA(fit,'shared',features=1:4)
local <- reconstruct_MSGPCA(fit,'local',features=1:4)
all <- reconstruct_MSGPCA(fit,features=1:4)
stopifnot(isTRUE(all.equal(all[[1]],shared[[1]]+local[[1]])))
scores <- get_scores(fit,'shared')
loadings <- get_loadings(fit,'shared')
head(sort(loadings[,1],decreasing=TRUE))
head(sort(loadings[,1],decreasing=FALSE))
head(sort(abs(loadings[,1]),decreasing=TRUE))
plot_factors(fit,'shared')
plot_factors(fit,'local')
plot_loadings(fit,top=8)

