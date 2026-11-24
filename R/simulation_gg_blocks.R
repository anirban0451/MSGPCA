#' Simulate one multi-slice ggblocks count dataset
#' @param n_slices Number of slices; centroid_shift requires at least five.
#' @param nside Positive grid side divisible by six.
#' @param Jsp,Jmix,Jns Numbers of purely spatial, mixed and nonspatial features.
#' @param Lns Number of Bernoulli nonspatial factors.
#' @param bkg_mean Nonnegative baseline mean.
#' @param nb_shape Positive negative-binomial size; variance is mu + mu^2/size.
#' @param expr_mean Nonnegative loading intensity.
#' @param mix_frac_spat Spatial share of mixed-gene loading intensity.
#' @param prob_nonspatial Bernoulli factor activation probability.
#' @param scenario One of 'none', 'slab', 'wave', 'centroid_shift'.
#' @param local_intensity Override scenario defaults: slab 20, wave 8, centroid 25.
#' @param local_slice Slice receiving slab or wave replacement.
#' @param mix_frac_red,mix_frac_green Fractions of background genes in disjoint
#'   replacement groups, rounded as in the original generator.
#' @param seed Optional seed with caller RNG state restored.
#' @return List with Y_list, locs_list, true_shared_factors, true_local_patterns,
#'   true_cluster_labels, generating loadings/factors, gene groups, perturbation
#'   masks and parameters. Baseline means are reconstructible from bkg_mean,
#'   true_shared_factors, true_spatial_loadings and nonspatial factors/loadings;
#'   perturbations specify changes to these means.
#' @details Preserves the original row-major block expansion and coordinate
#'   order (x increases first, y decreases by rows). No automatic plotting or
#'   files. Wave/slab replace counts in background genes; centroid shifts add
#'   independent counts to purely spatial genes in slices 3--5. Extra slices are
#'   unperturbed. Baseline cluster labels always describe the conserved domains.
#'   Defaults total 400 genes, unlike the manuscript's reported 500. Generating
#'   patterns/loadings are not orthogonal model parameters. Lns=0 is allowed only
#'   if Jns=0 and mixed loadings have spatial fraction one (or Jmix=0).
#' @examples
#' sim <- simulation_gg_blocks(nside=6,Jsp=8,Jmix=12,Jns=10,seed=101)
#' dim(sim$Y_list[[1]])
#' @export
simulation_gg_blocks <- function(n_slices=2L,nside=36L,Jsp=50L,Jmix=300L,Jns=50L,Lns=3L,
 bkg_mean=0.2,nb_shape=10,expr_mean=20,mix_frac_spat=0.55,prob_nonspatial=0.2,
 scenario='wave',local_intensity=NULL,local_slice=2L,mix_frac_red=0.2,mix_frac_green=0.2,seed=NULL) {
  scenario <- match.arg(scenario,c('none','slab','wave','centroid_shift'))
  for (nm in c('n_slices','nside')) .integer(get(nm),nm,1)
  for (nm in c('Jsp','Jmix','Jns','Lns')) .integer(get(nm),nm)
  if (nside%%6 || Jsp+Jmix+Jns==0) stop('nside must be divisible by six and total features positive')
  for (nm in c('mix_frac_spat','prob_nonspatial','mix_frac_red','mix_frac_green')) .scalar(get(nm),nm,0,1)
  if (mix_frac_red+mix_frac_green>1) stop('Perturbation fractions must sum to <=1')
  .scalar(bkg_mean,'bkg_mean'); .scalar(expr_mean,'expr_mean'); .scalar(nb_shape,'nb_shape')
  if(nb_shape==0) stop('nb_shape must be positive')
  if(!Lns && (Jns>0 || Jmix>0 && mix_frac_spat!=1)) stop('Zero nonspatial factors incompatible with feature composition')
  if (is.null(local_intensity)) local_intensity <- switch(scenario,none=0,slab=20,wave=8,centroid_shift=25)
  .scalar(local_intensity,'local_intensity')
  if(scenario %in% c('wave','slab')) {
    .integer(local_slice,'local_slice',1); if(local_slice>n_slices) stop('Invalid perturbation slice')
    if(Jns==0 && (mix_frac_red>0||mix_frac_green>0)) stop('Background genes required for perturbations')
  }
  if(scenario=='centroid_shift' && (n_slices<5||Jsp==0)) stop('Centroid shift requires >=5 slices and spatial genes')
  parameters <- as.list(environment()); parameters <- parameters[setdiff(names(parameters),c('nm'))]
  .with_seed(seed,function() {
    geo <- .gg_geometry(nside); F <- geo$F; z <- geo$coords; N <- nside^2; J <- Jsp+Jmix+Jns
    dummy <- function(l) {if (!l) return(matrix(0,J,0)); a <- sample(seq_len(l),J,replace=TRUE); diag(l)[a,,drop=FALSE]}
    W <- dummy(4); V <- dummy(Lns)
    sp <- seq_len(Jsp); mix <- Jsp+seq_len(Jmix); ns <- Jsp+Jmix+seq_len(Jns)
    W[sp,] <- W[sp,,drop=FALSE]*expr_mean; V[sp,] <- 0
    W[mix,] <- W[mix,,drop=FALSE]*mix_frac_spat*expr_mean
    V[mix,] <- V[mix,,drop=FALSE]*(1-mix_frac_spat)*expr_mean
    W[ns,] <- 0; V[ns,] <- V[ns,,drop=FALSE]*expr_mean
    rownames(W) <- rownames(V) <- paste0('Gene',seq_len(J))
    Y <- U <- coords <- vector('list',n_slices)
    for(m in seq_len(n_slices)) {
      U[[m]] <- matrix(stats::rbinom(N*Lns,1,prob_nonspatial),N,Lns)
      mu <- bkg_mean+tcrossprod(F,W)+tcrossprod(U[[m]],V)
      Y[[m]] <- matrix(stats::rnbinom(N*J,size=nb_shape,mu=mu),N,J,dimnames=list(rownames(z),rownames(W)))
      coords[[m]] <- z
    }
    masks <- list(); genes <- list(); affected <- integer(); mode <- 'none'
    if(scenario %in% c('wave','slab')) {
      red <- ns[seq_len(round(Jns*mix_frac_red))]
      green <- ns[round(Jns*mix_frac_red)+seq_len(round(Jns*mix_frac_green))]
      boundary <- if(scenario=='wave') 0.5+0.15*sin(4*pi*z[,1])+0.04*cos(10*pi*z[,1]) else rep(0.3,N)
      masks <- list(wave_or_slab=z[,2]<=boundary,diagonal=z[,1]+z[,2]<=0.9)
      genes <- list(wave_or_slab=red,diagonal=green); affected <- local_slice; mode <- 'replacement'
      for(a in seq_along(masks)) for(g in genes[[a]]) {
        idx <- which(masks[[a]])
        Y[[local_slice]][idx,g] <- stats::rnbinom(length(idx),size=nb_shape,mu=local_intensity)
      }
    }
    if(scenario=='centroid_shift') {
      masks <- list(wave=z[,2]<=0.5+0.15*sin(4*pi*z[,1])+0.04*cos(10*pi*z[,1]),
        diagonal=z[,1]+z[,2]>=1.2,circle=(z[,1]-0.5)^2+(z[,2]-0.5)^2<=0.08)
      genes <- rep(list(sp),3); names(genes) <- names(masks); affected <- 3:5; mode <- 'additive'
      for(a in seq_along(masks)) {
        idx <- which(masks[[a]]); m <- affected[a]
        Y[[m]][idx,sp] <- Y[[m]][idx,sp,drop=FALSE]+matrix(stats::rnbinom(length(idx)*length(sp),size=nb_shape,mu=local_intensity),length(idx),length(sp))
      }
    }
    labels <- rep('Background',N)
    for(l in seq_len(4)) labels[F[,l]==1] <- paste0('Domain_Factor_',l)
    sn <- paste0('Slice',seq_len(n_slices)); names(Y) <- names(U) <- names(coords) <- sn
    patterns <- lapply(seq_len(n_slices),function(m) {
      a <- matrix(0,N,length(masks),dimnames=list(rownames(z),names(masks)))
      if(mode=='replacement' && m==local_slice) for(j in seq_along(masks)) a[,j] <- masks[[j]]
      if(mode=='additive' && m %in% affected) a[,match(m,affected)] <- masks[[match(m,affected)]]
      a
    }); names(patterns) <- sn
    list(Y_list=Y,locs_list=coords,true_shared_factors=F,true_local_patterns=patterns,
      true_cluster_labels=stats::setNames(rep(list(labels),n_slices),sn),true_spatial_loadings=W,
      true_nonspatial_loadings=V,true_nonspatial_factors=U,gene_groups=list(spatial=sp,mixed=mix,nonspatial=ns),
      local_gene_indices=genes,perturbations=list(masks=masks,affected_slices=affected,intensity=local_intensity,mechanism=mode),parameters=parameters)
  })
}
