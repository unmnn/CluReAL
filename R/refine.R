#' Refine a clustering solution using the CluReAL algorithm
#'
#' Run the CluReAL V2 algorithm on a data matrix and a vector of cluster 
#' membership indices.
#' The output is a list with the refined cluster membership indices and the 
#' the cluster context. 
#' 
#' For a description of the algorithm, see:
#' 
#' Iglesias, Felix, Tanja Zseby, and Arthur Zimek. "Clustering refinement." 
#' International Journal of Data Science and Analytics (2021): 1-21.  
#' URL: \url{https://doi.org/10.1007/s41060-021-00275-z}.
#'
#' @param x The input data (numeric matrix).
#' @param y Cluster membership indices (integer vector). -1 indicates noise.
#' @param cc The \code{\link[cluster_context]{cluster context}}.
#' @param gv The \code{\link[gval]{GOI cluster validity}} indices object.
#' @param rc The \code{\link[refinement_context]{refinement context}}.
#' @param rep The number of refinement repetitions.
#' @param min_rdens The minimum density threshold (relative to the overall
#' density) of a cluster to not be considered noise.
#' @param min_mass The relative minimum size of a cluster to not be considered 
#' noise.
#' @param out_sens Outlier sensitivity. Specifies how far away outliers must be 
#' from the centroids to avoid trying to merge them with the clusters.
#' @param prun_level Pruning level. 0: normal (default), 
#' 1: high cluster overlap expected, 
#' 2: very high cluster overlap expected
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{y}: Vector of the refined cluster membership indices.
#'   \item \code{cc}: Refined cluster context
#' }
#' @export
#'
#' @examples
#' # TODO
refine <- function(
  x, y, cc, gv, rc, 
  rep = 0, 
  min_rdens = -0.8, 
  min_mass = 0.005, 
  out_sens = 0.75, 
  prun_level = 0
) {
  
  for(j in 1:(rep+1)) {
    if(any(rc$mm)) {
      y <- dig_multimodal(x, y, rc$mm)
      cc <- cluster_context(x, y)
      gv <- gval(cc)
      rc <- refinement_context(x, y, cc, gv)
    }
    
    ch_flag <- FALSE
    
    # transform objects of hazy or low-mass clusters into outliers
    for(i in 1:cc$k) {
      if(cc$outliers < 0.5 & 
         (rc$k_dens[i] <= min_rdens | cc$mass[i] < (sum(cc$mass) * min_mass))) {
        y[y == i] <- -1 # or reassign them
        ch_flag <- TRUE
      }
    }
    
    if(ch_flag) {
      y <- rebuilt_labels(y)
      cc <- cluster_context(x, y)
      gv <- gval(cc)
      rc <- refinement_context(x, y, cc, gv)
    }
    
    y <- graph_ref(x, y, rc$kinship, prun_level)
    
    cc <- cluster_context(x, y)
    gv <- gval(cc)
    rc <- refinement_context(x, y, cc, gv)
    
    if(any(y == -1)) {
      y <- reassign_outliers(x, y, out_sens, cc$centroids, gv$str_r)
      cc <- cluster_context(x, y)
    }
    
    if(j < rep) {
      gv <- gval(cc)
      rc <- refinement_context(x, y, cc, gv)
    }
    
  }
  
  list(y = y, cc = cc)
}
