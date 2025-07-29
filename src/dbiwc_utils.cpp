#include<RcppArmadillo.h>
#include <chrono>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(transport)]]
using namespace Rcpp;


/**
 * Compute the Wasserstein distance between two one-dimensional samples
 *
 * @param sam1 First sample vector (arma::vec)
 * @param sam2 Second sample vector (arma::vec)
 * @param p Power parameter for Wasserstein distance (default = 1)
 * @return The computed Wasserstein distance between sam1 and sam2
 */
static double wasserstein(const arma::vec& sam1, const arma::vec& sam2, int p = 1){
  // Get the wasserstein1d function from the transport package namespace
  Function wasserstein1d("wasserstein1d", Environment::namespace_env("transport"));

  // Compute the distance with specified power parameter p
  double dist = as<double>(wasserstein1d(sam1, sam2, Named("p") = p));

  return dist;
}


/**
 * Calculate cluster centers from grouped samples
 *
 * @param sam List of matrices, where each matrix represents a sample
 * @param k Integer, number of clusters
 * @param clus Vector of cluster assignments (1-based indexing)
 * @return List of matrices, where each matrix is the combined samples for a cluster center
 */
static List clus_center(List sam, int k, arma::vec clus) {
  // Initialize a list to store the k cluster centers
  List center(k);

  // Loop through each cluster
  for(int i = 0; i < k; i++) {
    // Find indices of samples belonging to current cluster (i+1 because clusters are 1-based)
    arma::vec clus_id = arma::conv_to<arma::vec>::from(find(clus == i + 1));
    int n_i = clus_id.size();  // Number of samples in current cluster

    // Calculate number of rows for each sample in this cluster
    arma::vec rows_i(n_i);
    for(int r = 0; r < n_i; r++) {
      rows_i(r) = as<arma::mat>(sam[clus_id(r)]).n_rows;
    }

    // Initialize matrix to store all samples from this cluster
    // Total rows = sum of rows from all samples in cluster
    // Columns same as input samples
    arma::mat center_i(sum(rows_i), as<arma::mat>(sam[0]).n_cols);

    // Combine all samples from this cluster into one matrix
    int rows_t = 0;  // Tracks current row position in center_i
    for(int t = 0; t < n_i; t++) {
      // Copy sample matrix into the appropriate rows of center_i
      center_i.rows(rows_t, rows_t + rows_i(t) - 1) = as<arma::mat>(sam[clus_id(t)]);
      // Update row position for next sample
      rows_t = rows_t + rows_i(t);
    }

    // Store combined matrix for this cluster
    center[i] = center_i;
  }

  // Return list of cluster centers
  return center;
}


/**
 * Calculate intra-cluster Wasserstein distance for a specific variable
 *
 * @param sam List of matrices containing sample data
 * @param clus Vector of cluster assignments (1-based indexing)
 * @param center Matrix representing the cluster center
 * @param g Integer index of the variable to analyze (0-based)
 * @param k0 Cluster label to analyze (1-based)
 * @return Average Wasserstein distance for the specified variable within the cluster
 */
static double dis_in(List sam, arma::vec clus, arma::mat center, int g, int k0) {
  // Get indices of all samples belonging to cluster k0
  arma::vec clus_id = arma::conv_to<arma::vec>::from(find(clus == k0));
  int ns = clus_id.size();  // Number of samples in cluster k0

  double sigma = 0;  // Accumulator for total intra-cluster distance

  // Calculate distance for each sample in the cluster
  for(int i = 0; i < ns; i++) {
    // Extract the sample matrix for current cluster member
    arma::mat sam_k0 = sam[clus_id[i]];

    // Compute Wasserstein distance between:
    // - The g-th variable of current sample
    // - The g-th variable of cluster center
    // The '1' parameter indicates we're using L1 distance (Earth Mover's Distance)
    double dis_i = wasserstein(sam_k0.col(g), center.col(g), 1);

    sigma = sigma + dis_i;  // Add to total distance
  }

  // Return average intra-cluster distance
  return sigma / ns;
}


/**
 * Calculate the Davies-Bouldin Index (DBI) for a specific variable
 *
 * DBI measures cluster separation quality, where lower values indicate better clustering.
 * This version computes DBI for a single variable (column) across all clusters.
 *
 * @param sam List of matrices, where each matrix represents a sample observation
 * @param clus Vector of cluster assignments (1-based indexing)
 * @param k Integer number of clusters
 * @param clus_center List of matrices representing cluster centers
 * @param g Integer index of the variable/column to analyze (0-based)
 * @return DBI value for the specified variable (lower values indicate better clustering)
 */
static double dbi_g(List sam, arma::vec clus, int k, List clus_center, int g) {
  // Create list to store indices of samples for each cluster
  List clus_id(k);
  for(int i = 0; i < k; i++) {
    clus_id[i] = find(clus == i + 1);  // Get indices for cluster i+1 (1-based)
  }

  // Initialize DBI matrix to store pairwise cluster comparisons
  arma::mat dbim(k, k);
  dbim.fill(0.0);  // Initialize all elements to 0

  // Calculate pairwise DBI between clusters
  for(int i = 0; i < k; i++) {
    for(int j = i + 1; j < k; j++) {
      // Compute inter-cluster distance (theta) using Wasserstein metric
      double theta = wasserstein(
        (as<arma::mat>(clus_center[i])).col(g),  // g-th variable of i-th cluster center
        (as<arma::mat>(clus_center[j])).col(g),  // g-th variable of j-th cluster center
        1  // Use L1 distance (Earth Mover's Distance)
      );

      // Compute intra-cluster distances (sigma) for both clusters
      double sigma_i = dis_in(sam, clus, clus_center[j], g, i + 1);
      double sigma_j = dis_in(sam, clus, clus_center[j], g, j + 1);

      // Store symmetric DBI ratio in matrix
      dbim(i, j) = theta / (sigma_i + sigma_j);
      dbim(j, i) = dbim(i, j);  // Matrix is symmetric
    }
  }

  // Find maximum DBI value for each cluster (worst-case separation)
  arma::vec dbimax = max(dbim, 1);  // Maximum of each row

  // Calculate final DBI as average of worst-case separations
  double dbi = sum(dbimax) / k;

  return dbi;
}


/**
 * Calculate weighted mixed-type distribution distance between two samples
 *
 * @param sam1 First sample matrix (rows: observations, columns: variables)
 * @param sam2 Second sample matrix (same dimensions as sam1)
 * @param w Weight vector for variables (length must match number of columns)
 * @return Weighted composite distance between the two samples
 */
static double dis_mixed(arma::mat sam1, arma::mat sam2, arma::vec w) {
  // Validate input dimensions
  if(sam1.n_cols != sam2.n_cols || sam1.n_cols != w.n_elem) {
    stop("Input dimension mismatch: samples and weights must have consistent dimensions");
  }

  int d = sam1.n_cols;       // Number of variables/dimensions
  double dis_mixed = 0;      // Initialize composite distance

  // Calculate weighted distance for each variable
  for(int i = 0; i < d; i++) {
    // Compute Wasserstein distance (Earth Mover's Distance) for current variable
    double dis_i = wasserstein(sam1.col(i), sam2.col(i), 1);  // L1 distance

    // Add weighted component to composite distance
    dis_mixed += dis_i * w(i);
  }

  return dis_mixed;
}


/**
 * Initialize clustering centers
 *
 * @param sam List of matrices where each matrix represents a sample observation
 * @param k Integer number of clusters to initialize
 * @return List of matrices representing initial cluster centers
 */
static List initialization(List sam, int k) {
  // Get data dimensions
  int d = (as<arma::mat>(sam[0])).n_cols;  // Number of variables/features
  int n = sam.size();                       // Total number of samples

  // Initialize uniform weights for distance calculation
  arma::vec w(d);
  w.fill(1/static_cast<double>(d));  // Equal weights for all variables

  // Create index vector of all sample IDs [0, 1, ..., n-1]
  arma::vec id = arma::regspace(0, 1, n-1);

  // Initialize container for cluster centers
  List center(k);

  // Step 1: Randomly select first cluster center
  int id0 = arma::randi<int>(arma::distr_param(0, n-1));
  center[0] = sam[id0];
  id = id.elem(find(id != id0));  // Remove selected sample from candidates

  // Step 2: Iteratively select remaining centers using weighted probabilities
  for(int i = 1; i < k; i++) {
    // Calculate distance matrix between candidates and existing centers
    arma::mat p(id.size(), i);  // Rows: candidate points, Columns: existing centers

    for(int r = 0; r < static_cast<int>(id.size()); r++) {
      for(int s = 0; s < i; s++) {
        // Compute mixed distance between candidate and existing center
        p(r,s) = dis_mixed(as<arma::mat>(sam[id(r)]),
          as<arma::mat>(center[s]),
          w);
      }
    }

    // Get minimum distance to nearest center for each candidate
    arma::vec pm = min(p, 1);  // Minimum along rows (nearest center)
    pm = pm / sum(pm);         // Convert to probabilities

    // Select new center weighted by squared distances
    int id_i = sample(as<NumericVector>(wrap(id)),  // Candidate indices
                      1,                           // Select 1 center
                      false,                       // Without replacement
                      as<NumericVector>(wrap(pm)))[0];  // Probabilities

    center[i] = sam[id_i];      // Store new center
    id = id.elem(find(id != id_i));  // Remove selected sample
  }

  return center;
}


/**
 * Perform DBI-WC Clustering Algorithm for mixed-type distribution data
 *
 * @param sam List of matrices where each matrix represents a sample observation
 * @param k Integer number of clusters
 * @param max_itr Maximum number of iterations (default: 20)
 * @param tol Convergence tolerance (default: 1e-4)
 * @return List containing:
 *         - label: Cluster assignments
 *         - mdbi: Davies-Bouldin Index values for each feature
 *         - weight: Learned feature weights
 *         - center: Final cluster centers
 *         - num_itr: Number of iterations performed
 */
// [[Rcpp::export]]
List DBIWC(List sam, int k, int max_itr = 20, double tol = 1e-4) {
  // Get data dimensions
  int d = (as<arma::mat>(sam[0])).n_cols;  // Number of features/variables
  int n = sam.size();                       // Total number of samples

  // Step 1: Initialize cluster centers using K-means++
  List center = initialization(sam, k);
  List center_old = center;  // Store for convergence checking

  // Initialize algorithm variables
  double loss = 1e6;         // Initialize with large loss value
  arma::vec loss_k(k);       // Per-cluster loss values
  int itr = 0;               // Iteration counter
  arma::mat dis(n, k);       // Distance matrix (samples x centers)
  arma::vec dbi(d);          // Davies-Bouldin Index values per feature
  arma::vec w(d);            // Feature weights
  w.fill(1/static_cast<double>(d));  // Initialize with uniform weights
  arma::vec clus(n);         // Cluster assignments (1-based)

  // Main algorithm loop
  while(loss > tol && itr < max_itr) {
    // Step 2: Calculate distances between samples and centers
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++) {
        dis(i,j) = dis_mixed(as<arma::mat>(sam[i]),
            as<arma::mat>(center[j]),
            w);
      }
    }

    // Step 3: Assign samples to nearest clusters
    clus = arma::conv_to<arma::vec>::from(index_min(dis, 1)) + 1;

    // Step 4: Update cluster centers
    center = clus_center(sam, k, clus);

    // Step 5: Update feature weights based on cluster quality
    for(int t = 0; t < d; t++) {
      dbi(t) = dbi_g(sam, clus, k, center, t);  // Calculate DBI per feature
    }
    w = dbi / sum(dbi);  // Normalize to get weights

    // Step 6: Check convergence
    for(int i = 0; i < k; i++) {
      loss_k(i) = dis_mixed(as<arma::mat>(center[i]),
             as<arma::mat>(center_old[i]),
             w);
    }

    // Relative change stopping criterion
    if(abs(loss - max(loss_k))/loss < tol) {
      itr = itr + 1;
      break;
    }

    loss = max(loss_k);     // Update overall loss
    center_old = center;    // Store current centers
    itr = itr + 1;          // Increment iteration counter
  }

  // Prepare results
  List res;
  res["label"] = clus;      // Cluster assignments (1-based)
  res["mdbi"] = dbi;        // Final DBI values per feature
  res["weight"] = w;        // Learned feature weights
  res["center"] = center;   // Final cluster centers
  res["num_itr"] = itr;     // Number of iterations performed

  return res;
}


/**
 * Performs SA-WC Clustering Algorithm for mixed-type distribution data
 *
 * @param sam List of matrices where each matrix represents a sample observation
 * @param k Number of clusters to identify (must be ≥ 2)
 * @param max_itr Maximum number of iterations allowed (default: 20)
 * @param tol Convergence tolerance threshold (default: 1e-4)
 * @return List containing:
 *         - label: Vector of cluster assignments (1-based)
 *         - weight: Vector of feature weights used
 *         - center: List of final cluster centers
 *         - num_itr: Number of iterations performed
 */
// [[Rcpp::export]]
List SAWC(List sam, int k, int max_itr = 20, double tol = 1e-4) {
  // Validate input dimensions
  if(k < 2) stop("Number of clusters k must be at least 2");

  // Get data dimensions
  int d = (as<arma::mat>(sam[0])).n_cols;  // Number of features/variables
  int n = sam.size();                      // Total number of samples

  // Step 1: Initialize cluster centers using K-means++
  List center = initialization(sam, k);
  List center_old = center;  // Store previous centers for convergence check

  // Initialize algorithm variables
  double loss = 1e6;         // Initialize with large loss value
  arma::vec loss_k(k);       // Per-cluster loss values
  int itr = 0;               // Iteration counter
  arma::mat dis(n, k);       // Distance matrix (n samples × k centers)
  arma::vec w(d);            // Fixed feature weights
  w.fill(1/static_cast<double>(d));  // Uniform weights (1/d for each feature)
  arma::vec clus(n);         // Cluster assignments (1-based indexing)

  // Main clustering loop
  while(loss > tol && itr < max_itr) {
    // Step 2: Calculate distances between all samples and centers
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < k; j++) {
        // Compute weighted mixed distance using fixed weights
        dis(i,j) = dis_mixed(as<arma::mat>(sam[i]),
            as<arma::mat>(center[j]),
            w);
      }
    }

    // Step 3: Assign each sample to nearest cluster
    clus = arma::conv_to<arma::vec>::from(index_min(dis, 1)) + 1;

    // Step 4: Update cluster centers
    center = clus_center(sam, k, clus);

    // Step 5: Check convergence by measuring center movement
    for(int i = 0; i < k; i++) {
      loss_k(i) = dis_mixed(as<arma::mat>(center[i]),
             as<arma::mat>(center_old[i]),
             w);
    }

    // Relative change stopping criterion
    if(abs(loss - max(loss_k))/loss < tol) {
      itr++;  // Count final iteration
      break;  // Exit if centers stabilized
    }

    // Update tracking variables
    loss = max(loss_k);     // Update overall loss
    center_old = center;    // Store current centers
    itr++;                  // Increment iteration counter
  }

  // Prepare output results
  List res;
  res["label"] = clus;      // Final cluster assignments (1-based)
  res["weight"] = w;        // Feature weights used (fixed)
  res["center"] = center;   // Final cluster centers
  res["num_itr"] = itr;     // Total iterations performed

  return res;
}
