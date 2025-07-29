#' DBI-WC Clustering Algorithm
#'
#' @title DBI-WC Clustering Algorithm
#' @description The DBI-WC clustering algorithm is a novel method designed for mixed-type distribution data. It employs a modified Davies-Bouldin Index (DBI) to construct a DBI-weighted Wasserstein distance metric, extending the traditional Wasserstein distance. Combined with a robust initialization strategy, DBI-WC achieves strong performance on mixed-type distribution data.
#'
#' @param sam List of matrices where each matrix represents a sample observation
#' @param k Integer number of clusters
#' @param max_itr Maximum number of iterations (default: 20)
#' @param tol Convergence tolerance (default: 1e-4)
#'
#' @return List containing:
#' \tabular{ll}{
#' \code{"label"} \tab Cluster assignments. \cr
#' \code{"mdbi"} \tab Davies-Bouldin Index values for each feature. \cr
#' \code{"weight"} \tab Learned distance weights. \cr
#' \code{"center"} \tab Final cluster centers. \cr
#' \code{"num_itr"} \tab Number of iterations performed. \cr
#' }
#' @export
dbiWC <- function(sam, k, max_itr=20, tol=1e-4){
  res <- DBIWC(sam, k, max_itr, tol)
  return(res)
}


#' @title SA-WC Clustering Algorithm
#' @description The SA-WC algorithm is a simplified version of DBI-WC. It uses a straightforward weighting scheme to integrate various distribution data. All other specifications are consistent with those of DBI-WC.
#'
#' @param sam List of matrices where each matrix represents a sample observation
#' @param k Integer number of clusters
#' @param max_itr Maximum number of iterations (default: 20)
#' @param tol Convergence tolerance (default: 1e-4)
#'
#' @return List containing:
#' \tabular{ll}{
#' \code{"label"} \tab Cluster assignments. \cr
#' \code{"weight"} \tab Learned distance weights. \cr
#' \code{"center"} \tab Final cluster centers. \cr
#' \code{"num_itr"} \tab Number of iterations performed. \cr
#' }
#' @export
saWC <- function(sam, k, max_itr=20, tol=1e-4){
  res <- SAWC(sam, k, max_itr, tol)
  return(res)
}
