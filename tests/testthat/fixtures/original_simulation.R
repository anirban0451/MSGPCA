# Verbatim supplied generator functions, retained for equivalence testing.
gen_ggblocks_base <- function() {
  A <- matrix(0, nrow = 4, ncol = 36)
  A[1, c(1, 6, 7, 8, 13) + 1] <- 1
  A[2, c(3, 4, 5, 9, 11, 15, 16, 17) + 1] <- 1
  A[3, c(18, 24, 25, 30, 31, 32) + 1] <- 1
  A[4, c(21, 22, 23, 28, 34) + 1] <- 1
  return(A)
}

gen_spatial_factors <- function(nside = 36) {
  A_base <- gen_ggblocks_base()
  unit <- 6
  if (nside %% unit != 0) stop(sprintf("nside must be divisible by %d", unit))
  ncopy <- nside / unit
  
  L <- nrow(A_base)
  N <- nside^2
  
  F_mat <- matrix(0, nrow = N, ncol = L)
  for (l in 1:L) {
    block_2d <- matrix(A_base[l, ], nrow = unit, ncol = unit, byrow = TRUE)
    inflated_2d <- block_2d %x% matrix(1, nrow = ncopy, ncol = ncopy)
    F_mat[, l] <- as.vector(t(inflated_2d))
  }
  return(F_mat)
}

gen_spatial_coords <- function(nside) {
  grid <- expand.grid(x = 1:nside, y = 1:nside)
  grid$y <- -grid$y # Flip Y-axis to mirror spatial plotting orientation
  
  # Standardize scales cleanly to [0, 1]
  grid$x <- (grid$x - min(grid$x)) / (max(grid$x) - min(grid$x))
  grid$y <- (grid$y - min(grid$y)) / (max(grid$y) - min(grid$y)) # <-- Fixed typo here
  
  return(as.matrix(grid))
}

gen_loadings <- function(Lsp, Lns = 3, Jsp = 0, Jmix = 500, Jns = 0, 
                         expr_mean = 20.0, mix_frac_spat = 0.55) {
  J <- Jsp + Jmix + Jns
  get_dummies <- function(vec, levels_count) {
    fac <- factor(vec, levels = 1:levels_count)
    return(model.matrix(~ fac - 1))
  }
  
  W <- if (Lsp > 0) get_dummies(sample(1:Lsp, J, replace = TRUE), Lsp) else matrix(0, nrow = J, ncol = 0)
  V <- if (Lns > 0) get_dummies(sample(1:Lns, J, replace = TRUE), Lns) else matrix(0, nrow = J, ncol = 0)
  
  if (Jsp > 0) {
    W[1:Jsp, ] <- W[1:Jsp, ] * expr_mean
    V[1:Jsp, ] <- 0
  }
  if (Jmix > 0) {
    mix_idx <- (Jsp + 1):(Jsp + Jmix)
    W[mix_idx, ] <- W[mix_idx, ] * (mix_frac_spat * expr_mean)
    V[mix_idx, ] <- V[mix_idx, ] * ((1 - mix_frac_spat) * expr_mean)
  }
  if (Jns > 0) {
    ns_idx <- (Jsp + Jmix + 1):J
    W[ns_idx, ] <- 0
    V[ns_idx, ] <- V[ns_idx, ] * expr_mean
  }
  
  return(list(W = W, V = V))
}

simulate_dataset <- function(n_slices = 2, nside = 36, Jsp = 50, Jmix = 300, 
                             Jns = 50, Lns = 3, bkg_mean = 0.2, nb_shape = 10.0, 
                             expr_mean = 20.0, mix_frac_spat = 0.55) {
  
  N <- nside^2
  J <- Jsp + Jmix + Jns
  
  F_spatial <- gen_spatial_factors(nside = nside)
  Lsp <- ncol(F_spatial)
  coords <- gen_spatial_coords(nside)
  
  shared_loadings <- gen_loadings(Lsp = Lsp, Lns = Lns, Jsp = Jsp, Jmix = Jmix, Jns = Jns,
                                  expr_mean = expr_mean, mix_frac_spat = mix_frac_spat)
  
  Ylist_out    <- list()
  locslist_out <- list()
  
  for(s in 1:n_slices) {
    U_nonspatial <- matrix(rbinom(N * Lns, size = 1, prob = 0.2), nrow = N, ncol = Lns)
    Lambda <- bkg_mean + (F_spatial %*% t(shared_loadings$W)) + (U_nonspatial %*% t(shared_loadings$V))
    Y_counts <- matrix(rnbinom(N * J, size = nb_shape, mu = Lambda), nrow = N, ncol = J)
    
    colnames(Y_counts) <- paste0("Gene", 1:J)
    rownames(Y_counts) <- paste0("Spot", 1:N)
    colnames(coords)   <- c("x", "y")
    rownames(coords)   <- paste0("Spot", 1:N)
    
    Ylist_out[[s]]    <- Y_counts
    locslist_out[[s]] <- coords
  }
  
  return(list(Ylist = Ylist_out, locslist = locslist_out, true_spfac = F_spatial))
}

inject_targeted_local_wave_features <- function(Ylist, locslist, F_spatial, slice_id = 2, total_bg_genes = 50,
                                           mix_frac_red = 0.20, mix_frac_green = 0.20,
                                           intensity_mean = 20.0, nb_shape = 10.0) {
  
  Y_mat <- Ylist[[slice_id]]
  coords <- locslist[[slice_id]]
  J <- ncol(Y_mat)
  N <- nrow(Y_mat)
  
  bg_gene_start <- J - total_bg_genes + 1
  bg_gene_end   <- J
  bg_gene_indices <- bg_gene_start:bg_gene_end
  
  n_red   <- round(total_bg_genes * mix_frac_red)
  n_green <- round(total_bg_genes * mix_frac_green)
  
  red_genes   <- bg_gene_indices[1:n_red]
  green_genes <- bg_gene_indices[(n_red + 1):(n_red + n_green)]
  
  # ----------------------------------------------------------------------------
  # THE FIX: REPLACING THE SLAB WITH OPTION 2 WAVY GEOMETRY
  # ----------------------------------------------------------------------------
  # Coords are in [0, 1]. We map x to a 2*pi scale to ensure proper wave frequency.
  # Adjusting terms to match your undulating sketch:
  # - Baseline vertical height: 0.5
  # - Low-frequency wave amplitude: 0.15 (controls macro shape)
  # - High-frequency ripple amplitude: 0.04 (controls smaller organic bumps)
  x_scaled <- coords[, "x"] * 2 * pi
  wave_boundary <- 0.5 + (0.15 * sin(x_scaled * 2)) + (0.04 * cos(x_scaled * 5))
  
  # Spots situated above this non-linear boundary profile (adjusting for y-axis orientation)
  red_spots <- which(coords[, "y"] <= wave_boundary)
  
  # Retain original diagonal triangle features
  green_spots <- which((coords[, "y"] + coords[, "x"]) <= 0.9)
  
  if(length(red_spots) > 0 && length(red_genes) > 0) {
    for(g in red_genes) {
      Y_mat[red_spots, g] <- rnbinom(length(red_spots), size = nb_shape, mu = intensity_mean)
    }
  }
  
  if(length(green_spots) > 0 && length(green_spots) > 0) { # Adjusted variable check safety
    for(g in green_genes) {
      Y_mat[green_spots, g] <- rnbinom(length(green_spots), size = nb_shape, mu = intensity_mean)
    }
  }
  
  Ylist[[slice_id]] <- Y_mat
  
  # ----------------------------------------------------------------------------
  # GENERATE TRUE LABELS (Strictly 5 Biological Domains For Both Slices)
  # ----------------------------------------------------------------------------
  bio_clusters <- rep("Background", N)
  for (l in 1:ncol(F_spatial)) {
    bio_clusters[F_spatial[, l] == 1] <- paste0("Domain_Factor_", l)
  }
  
  return(list(
    Ylist = Ylist,
    true_cluster_labels = list(bio_clusters, bio_clusters)
  ))
}
