#' Create the cluster context
#'
#' @inheritParams refine
#'
#' @return A list with the following elements: 
#' \itemize{
#'   \item \code{k}: The number of clusters (integer).
#'   \item \code{centroids}: Centroid coordinates (k x n matrix).
#'   \item \code{mass}: Cluster sizes (integer vector of length k).
#'   \item \code{mn_da}: Mean Euclidean distances of the cluster members to 
#'   their centroid (double vector of length k).
#'   \item \code{md_da}: Median Euclidean distances of the cluster members to 
#'   their centroid (double vector of length k).
#'   \item \code{sd_da}: Standard deviation of the Euclidean distances of the 
#'   cluster members to their centroid (double vector of length k).
#'   \item \code{de}: Euclidean distance matrix for the centroids (k x k 
#'   matrix).
#'   \item \code{outliers}: The number of outliers (integer).
#' }
#' @export
#'
#' @examples
#' # TODO
cluster_context <- function(x, y) {

  k <- max(y)
  cc <- list(
    k = k,
    # centroids: matrix (kxm) with robust cluster centers
    centroids = matrix(0, nrow = k, ncol = ncol(x)), 
    # mass: cluster mass of cardinality (vector of length k)
    mass = vector("double", length = k), 
    # mn_da: cluster mean intra distance (vector of length k)
    mn_da = vector("double", length = k), 
    # md_da: cluster median intra distance (vector of length k)
    md_da = vector("double", length = k),
    # sd_da: cluster std-dev intra distance (vector of length k)
    sd_da = vector("double", length = k),
    # de: cluster inter distance matrix (kxk) 
    de = matrix(0, nrow = k, ncol = k), 
    # outliers: number of outliers / total data points (double)
    outliers = 0
  )
  
  c_x_i <- matrix(0, nrow = k, ncol = ncol(x))
  
  for(i in 1:k) {
    x_i <- x[y == i, , drop = FALSE]
    cc$mass[i] <- nrow(x_i)
    c_x <- apply(x_i, 2, median)
    c_x_i[i, ] <- c_x
    # dm: distances of cluster members to centroid
    dm <- unname(as.matrix(dist(rbind(c_x, x_i)))[1, ]) 
    cc$mn_da[i] <- mean(dm)
    cc$md_da[i] <- median(dm)
    cc$sd_da[i] <- sd(dm)
  }
  cc$de <- as.matrix(dist(c_x_i))
  cc$centroids <- c_x_i
  cc$outliers <- sum(y == -1) / sum(cc$mass)
  
  cc
} 