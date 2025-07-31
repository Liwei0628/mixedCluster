# mixedCluster


## Description
The R package `mixedCluster` implements two clustering algorithms for mixed-type distribution data, referred to as ***DBI-WC*** and ***SA-WC***, respectively. The **DBI-WC** method employs a modified Davies-Bouldin Index (DBI) to develop a DBI-weighted Wasserstein distance metric, generalizing the traditional Wasserstein distance to accommodate mixed-type distribution data. Coupled with a robust initialization strategy, ***DBI-WC*** typically achieves strong performance on mixed-type distribution data. ***SA-WC*** offers a simplified alternative with a straightforward weighting scheme, while maintaining the overall clustering framework.


## Installation
You can install this R package by running the following code in your `R` console:

```R
# R package `transport` is needed
install.packages("transport")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("Liwei0628/mixedCluster")
```


## Examples
Next, we will present some concrete examples to assist you in utilizing two main functions, `saWC` and `dbiWC`, from the R package `mixedCluster`, where `saWC` refers to the ***SA-WC*** method and `dbiWC` corresponds to the ***DBI-WC*** method.

#### 0. Load required R packages
You should first load the following R packages. If any are missing, please install them.

```R
library(mixedCluster)   # Our R package
library(aricode)        # clustering evaluation
```

#### 1. Generate simulated data
The simulated generation of mixed-type distribution data is implemented using the following R code:

```R
set.seed(666)
K = 3  # Number of cluster centers
N0 = 10  # Number of samples per cluster
lam = 25  # Constant related to number of observations per sample

# Function to generate simulated data
gnrt_data <- function(K, N0, lam){
  sam_ls <- list()  # Initialize list to store samples
  for(i in 1:K){
     for(j in 1:N0){
       # Generate number of observations for each sample
       vij = floor(max(rnorm(n = 1, mean = lam, sd = 0.1*lam), log(K*N0)))
       # Generate sample data based on cluster
       Xij <- switch(as.character(i),
                     # Cluster 1: Normal + Beta + Binomial
                      "1" = c(rnorm(vij, 5, 1), rbeta(vij, 3, 2), rbinom(vij, 1, 2/7)),
                     # Cluster 2: Normal + Gamma + Binomial
                      "2" = c(rnorm(vij, 2, 1), rgamma(vij, 3, 1), rbinom(vij, 1, 1/2)),
                     # Cluster 3: Normal + Lognormal + Binomial
                      "3" = c(rnorm(vij, 7, 1), rlnorm(vij, 5, 2), rbinom(vij, 1, 1/2))
      )
       Xij <- matrix(Xij, nrow = vij)  # Convert to matrix
       sam_ls[[j+N0*(i-1)]] <- Xij  # Store sample
     }
  }
  labels <- rep(1:4, each = N0)  # Create cluster labels
  clus_ls <- list("sam"=sam_ls, "label"=labels)  # Combine samples and labels
  return(clus_ls)
}

# Generate clustered data
clus_data = gnrt_data(K, N0, lam)
```

#### 2. SA-WC clustering algorithm
We can implement the ***SA-WC*** method for simulated mixed-type distribution data using the `saWC` function, as follows

```R
# Extract samples for clustering from the generated data
sam = clus_data$sam  # Samples to be clustered

# Basic Settings
k_pre = 4            # Predefined number of clusters
max_itr = 20         # Maximum number of iterations
tol = 1e-4           # Convergence tolerance threshold

# Run SA-WC clustering algorithm
res_SAWC <- saWC(sam, k_pre, max_itr, tol)

# View results:
res_SAWC$label       # Cluster assignments for each sample
res_SAWC$weight      # Distance weights for each variable (equal weights)
# res_SAWC$center    # Final cluster centers
res_SAWC$num_itr     # Number of iterations until convergence
```

#### 3. DBI-WC clustering algorithm
We can implement the ***DBI-WC*** methodfor simulated mixed-type distribution data using the `dbiWC` function, as follows

```R
# Extract samples for clustering from the generated data
sam = clus_data$sam  # Samples to be clustered

# Basic Settings
k_pre = 4            # Predefined number of clusters
max_itr = 20         # Maximum number of iterations
tol = 1e-4           # Convergence tolerance threshold

# Run DBI-WC clustering algorithm
res_DBIWC <- dbiWC(sam, k_pre, max_itr, tol)

# View results:
res_DBIWC$label       # Cluster assignments for each sample
res_DBIWC$mdbi        # Modified DBI values for each feature
res_DBIWC$weight      # Distance weights for each variable
# res_DBIWC$center    # Final cluster centers
res_DBIWC$num_itr     # Number of iterations until convergence
```


## References
Lin, L., Chen, Y., Li, H., Li, Y., & Wang, F. <sup>*</sup> (2025+). Clustering Taxi Service Patterns Using a DBI-Weighted Wasserstein Distance for Mixed-Type Distribution Data. Submitted.
