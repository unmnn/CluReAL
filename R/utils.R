rdensity <- function(x, y, k) {
  # inputs
  #   X: dataset (nxm), matrix of n vectors with m dimensions
  #   y: array with cluster labels (-1 for outliers) 
  #   k: number of clusters
  # outputs
  #   global_c_dens: global/overall density (scalar)
  #   k_dens: cluster relative densities (kx1-array)
  
  global_c <- apply(x, 2, median)
  
  d_x_to_global_c <- unname(as.matrix(dist(rbind(global_c, x))))[1, ]
  
  global_c_dens <- 1 / ((mean(d_x_to_global_c) + 2 * sd(d_x_to_global_c)) + nrow(x))
  
  k_dens <- vector("double", k)
  
  for(i in 1:k) {
    x_i <- x[y == i, , drop = FALSE]
    c_x_i <- apply(x_i, 2, median)
    intra_x_i <- unname(as.matrix(dist(rbind(c_x_i, x_i)))[1, ])
    medin_x_i <- median(intra_x_i)
    if(medin_x_i == 0) medin_x_i <- 1 # !! Why??
    i_card <- sum(y == i)
    k_dens[i] <- -1 + (i_card/medin_x_i)/global_c_dens
  }
  
  return(list(global_c_dens = global_c_dens, k_dens = k_dens))
}


dig_multimodal <- function(x, y, mm) {
  k <- max(y)
  c_id_offset <- 0
  n_clusters <- 2
  alg <- function(x) ClusterR::MiniBatchKmeans(x, clusters = n_clusters, seed = 10)
  for(i in 1:k) {
    if(mm[i]) {
      x_i <- x[y == i, , drop = FALSE]
      clustering <- alg(x_i)
      y_i <- ClusterR::predict_KMeans(x_i, CENTROIDS = clustering$centroids)
      d <- max(y_i) - 1
      y_i[y_i > 1] <- k + c_id_offset + d
      y_i[y_i == 1] <- i
      y[y == i] <- y_i
      c_id_offset <- c_id_offset + d
    }
  }
  y
}

cluster_kinship <- function(k, de, erad, srad) {
  # inputs
  #   k: number of clusters
  #   De: cluster inter distance matrix (k x k matrix)
  #   erad: extended radii (k-size array)
  #   srad: strict radii (k-size array)
  # outputs
  #   kinship: cluster kinship indices (k x k matrix): 
  #            5-unrelated, 
  #            4-acquaintances, 
  #            3-close-friends, 
  #            2-relatives, 
  #            1-parent and child, 
  #            0-itself.
  
  kinship <- matrix(0, nrow = k, ncol = k)
  # comb <- unname(as.matrix(expand.grid(1:k, 1:k)))
  comb <- t(combn(k, 2))
  for(r in 1:nrow(comb)) {
    # for(j in 1:k) {
    i <- comb[r, 1]
    j <- comb[r, 2]
    if(erad[i] + erad[j] <= de[i,j]) {
      # unrelated
      kinship[i,j] <- kinship[j,i] <- 5
    } else {
      if((erad[i] < de[i,j]) & (erad[j] < de[i,j])) {
        if(((de[i,j] - srad[i]) < erad[j]) | ((de[i,j] - srad[j]) < erad[i])) {
          # close friends
          kinship[i,j] <- kinship[j,i] <- 3
        } else {
          # acquaintance
          kinship[i,j] <- kinship[j,i] <- 4
        }
      } else if(((erad[i] + de[i,j]) < erad[j]) | 
                ((erad[j] + de[i,j]) < erad[i])){
        # parent and child
        kinship[i,j] <- kinship[j,i] <- 1
      } else {
        # relatives
        kinship[i,j] <- kinship[j,i] <- 2
      }
    }
    # }
  }
  kinship
}

multimodal_clusters <- function(x, y, k) {
  # inputs
  #   X: dataset (nxm), matrix of n vectors with m dimensions
  #   y: array with cluster labels (-1 for outliers) 
  #   k: number of clusters
  # outputs
  #   mm: multimodality flags for each cluster (kx1-array)
  mm <- vector("logical", length = k)
  for(i in 1:k) {
    x_i <- x[y == i, , drop = FALSE]
    if(nrow(x_i) > 0) {
      mm[i] <- multimodality(x_i)
    }
  }
  
  mm
}

multimodality <- function(x_i) {
  # inputs
  #   Xi: cluster data (nxm), matrix of n vectors with m dimensions
  # outputs
  #   mm: multimodality flag (scalar: 0 or 1)
  
  mm <- FALSE
  bwf <- 8
  
  for(i in 1:ncol(x_i)) {
    feat <- x_i[ , i]
    bw <- (max(feat) - min(feat)) / bwf
    if(bw > 0) {
      kde <- density(feat, adjust = bw) # !! leaf_size = 100???
      xbasis <- seq(min(feat), max(feat), length.out = 5 * nrow(x_i))
      x_kde <- approx(kde$x, kde$y, xbasis)$y
      peaks <- which(diff(sign(diff(x_kde))) == -2)
      
      if(length(peaks) > 1) {
        mm <- TRUE
        break
      }
    }
  }
  
  mm
}

rebuilt_labels <- function(y) {
  # inputs
  #   y: array with cluster labels (-1 for outliers) 
  # outputs
  #   y_new: refined array with cluster labels (-1 for outliers) 
  y_rem <- unique(y)
  outs <- which(y_rem == -1)
  a <- 0
  if(length(outs) > 0) {
    a <- 1
  }
  y_new <- y
  for(i in 1:(length(y_rem)-a)) {
    y_new[y == y_rem[i+a]] <- i
  }
  
  y_new
}

graph_ref <- function(x, y, kinship, prun_level) {
  kinship[kinship == 5] <- 0 # no edge between unrelated nodes
  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = kinship, mode = "undirected", weighted = TRUE
  )
  
  # kin <- NULL
  edgelist <- igraph::as_edgelist(g)
  # for(i in 1:nrow(edgelist)) {
  # kin <- c(kin, 5 - kinship[edgelist[i, 1], edgelist[i, 2]])
  # }
  
  # pos <- igraph::layout_with_fr(g)
  for(i in 1:nrow(edgelist)) {
    edge <- igraph::E(g)[[igraph::`%--%`(edgelist[i, 1], edgelist[i, 2])]]
    if(edge$weight >= (4 - prun_level)) {
      g <- igraph::delete.edges(g, edge)
    } else if(edge$weight == (3 - prun_level)) {
      if(multimodality(rbind(x[y == igraph::tail_of(g, edge), ], 
                             x[y == igraph::head_of(g, edge), ]))) {
        g <- igraph::delete.edges(g, edge)
      }
    }
  }
  
  lsub_g <- igraph::components(g)
  
  if(lsub_g$no == 1) {
    edgelist <- igraph::as_edgelist(g)
    for(i in 1:nrow(edgelist)) {
      edge <- igraph::E(g)[[igraph::`%--%`(edgelist[i, 1], edgelist[i, 2])]]
      if(edge$weight >= 2) {
        g <- igraph::delete.edges(g, edge) 
        # !! Why delete kinship=2 edges (relatives) only now and not 
        # already in previous for-loop?
      }
    }
  }
  
  lsub_g <- igraph::components(g)
  
  ynew <- vector("integer", length(y))
  ynew[y == -1] <- -1
  
  nc <- 1
  for(sub_g in 1:max(lsub_g$membership)){
    sub_g_idx <- which(lsub_g$membership == sub_g)
    for(lab in sub_g_idx) {
      ynew[y == lab] <- nc
    }
    nc <- nc + 1
  }
  
  ynew
}

reassign_outliers <- function(x, y, out_sens, centroids, ext_r) {
  if(out_sens == 0) {
    mem_th <- rep(Inf, length(ext_r))
  } else {
    mem_th <- ext_r / out_sens
  }
  
  x_i <- x[y == -1, , drop = FALSE]
  dm <- unname(as.matrix(dist(rbind(centroids, x_i)))[1:nrow(centroids), 
                                                      -(1:nrow(centroids)), 
                                                      drop = FALSE])
  y_out <- apply(dm, 2, which.min)
  dm_min <- apply(dm, 2, min)
  ths_ext <- mem_th[y_out]
  y_out[dm_min > ths_ext] <- -1
  y_new <- y
  y_new[y == -1] <- y_out
  y_new
  
}