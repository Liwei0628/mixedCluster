# mixedCluster


## Description
The R package `mixedCluster` implements two clustering algorithms for mixed-type distribution data, referred to as **DBI-WC** and **SA-WC**, respectively. The **DBI-WC** method employs a modified Davies-Bouldin Index (DBI) to develop a DBI-weighted Wasserstein distance metric, generalizing the traditional Wasserstein distance to accommodate mixed-type distribution data. Coupled with a robust initialization strategy, **DBI-WC** typically achieves strong performance on mixed-type distribution data. **SA-WC** offers a simplified alternative with a straightforward weighting scheme, while maintaining the overall clustering framework.


## Installation
You can install this R package by running the following code in your `R` console:

```R
# R package `transport` is needed
install.packages("transport")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("Liwei0628/mixedCluster")
```


## Examples
Next, we will present some concrete examples to assist you in utilizing two main functions, `saWC` and `dbiWC`, from the R package `mixedCluster`, where `saWC` refers to the **SA-WC** method and `dbiWC` corresponds to the **DBI-WC** method.
