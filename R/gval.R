#' Compute GOI cluster validity indices
#'
#' @inheritParams refine
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item \code{g_str}: Strict G-index (double).
#'   \item \code{g_rex}: Relaxed G-index (double).
#'   \item \code{g_min}: Min G-index (double)-
#'   \item \code{oi_st}: Cluster-individual strict overlap indices (double 
#'   vector of length k).
#'   \item \code{oi_rx}: Cluster-individual relaxed overlap indices (double 
#'   vector of length k).
#'   \item \code{oi_mn}: Cluster-individual min overlap indices (double vector 
#'   of length k).
#'   \item \code{ext_r}: Extended cluster radii (double vector of length k).
#'   \item \code{str_r}: Strict cluster radii (double vector of length k).
#'   \item \code{vol_r}: Extended-to-core ratio.
#' }
#' @export
#'
#' @examples
#' # TODO
gval <- function(cc) {

  k <- cc$k
  gv <- list(
    g_str = 0, # strict global index
    g_rex = 0, # relaxed global index
    g_min = 0, # minimum global index
    oi_st = vector("double", k), # indiv. strict indices
    oi_rx = vector("double", k), # indiv. relaxed indices 
    oi_mn = vector("double", k), # indiv. min indices
    ext_r = vector("double", k), # extended radii
    str_r = vector("double", k), # strict radii
    vol_r = vector("double", k) # times that the extended radius is in the core radius
  )
  
  radm <- matrix(0, nrow = k, ncol = 1)
  radm2 <- radm
  oi_st <- matrix(Inf, nrow = k, ncol = k)
  oi_rx <- oi_st
  
  gv$ext_r <- cc$mn_da + 2 * cc$sd_da
  gv$str_r <- cc$md_da
  gv$vol_r <- gv$ext_r / gv$str_r
  gv$vol_r[is.infinite(gv$vol_r)] <- 0
  
  radm <- gv$ext_r * cc$mass
  radm2 <- gv$str_r * cc$mass
  
  # per <- t(combn(1:k, m = 2))
  # per <- unname(as.matrix(expand.grid(1:k, 1:k)))!!!
  per <- t(combn(k, 2))
  per <- rbind(per, per[,c(2,1)])
  
  for(i in 1:nrow(per)) {
    # for(j in 1:ncol(per)) {
    n1 <- per[i,1]
    n2 <- per[i,2]
    oi_st[n1,n2] <- cc$de[n1,n2] - gv$ext_r[n1] - (cc$mn_da[n2] + 2 * cc$sd_da[n2]) # !! replace last term with gv$ext_r[n2] ???
    oi_rx[n1,n2] <- cc$de[n1,n2] - cc$md_da[n1] - cc$md_da[n2]
    # }
  }
  
  gv$oi_st <- apply(oi_st, 1, min)
  gv$oi_rx <- apply(oi_rx, 1, min)
  gv$oi_mn <- gv$oi_st / gv$ext_r
  gv$oi_mn[is.infinite(gv$oi_mn)] <- 0
  
  gv$g_str <- sum(gv$oi_st * cc$mass) / sum(radm)
  gv$g_rex <- sum(gv$oi_rx * cc$mass) / sum(radm2)
  
  if(length(gv$oi_mn)) {
    gv$g_min <- min(gv$oi_mn)
  }
  
  # if(verbose) {
  #   cat("- Validity index > GOI > Grex: ", gv$g_rex, "\n")
  #   cat("- Validity index > GOI > Gstr: ", gv$g_str, "\n")
  #   cat("- Validity index > GOI > Gmin: ", gv$g_min, "\n")
  # }
  
  gv
  
}