#' Create the refinement context
#'
#' @inheritParams refine 
#'
#' @return A list with four elements:
#' \itemize{
#'   \item \code{mm}: Cluster multimodality flag (logical vector of length k).
#'   \item \code{k_dens}: Cluster-individual relative density (double vector of 
#'   length k).
#'   \item \code{global_c_dens}: Global density (double).
#'   \item \code{kinship}: cluster kinship matrix (k x k matrix): 0-itself, 
#'   1-parent and child, 2-relatives, 3-close friends, 4-acquaintances, 
#'   5-unrelated. 
#' }
#' @export
#'
#' @examples
#' # TODO
refinement_context <- function(x, y, cc, gv) {
  rc <- list(
    mm = vector("logical", cc$k), # multimodality flags for each cluster
    k_dens = vector("double", cc$k), # cluster-relative densities
    global_c_dens = vector("double", cc$k), # global/overall density
    kinship = matrix(0, nrow = cc$k, ncol = cc$k) #cluster kinship indices
  )
  
  temp <- rdensity(x, y, cc$k)
  rc$k_dens <- temp$k_dens
  rc$global_c_dens <- temp$global_c_dens
  rc$kinship <- cluster_kinship(cc$k, cc$de, gv$ext_r, gv$str_r)
  
  rc
  
}