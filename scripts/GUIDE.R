
library(fastICA) ## package for implementation of fastICA


## Functions for the decomposition of a GWAS summary statistics matrix using the method
## Genetic Unmixing by Independent Decomposition (GUIDE).

## Function to get GUIDE unmixing matrices from multiple ICA runs

get.unmixing.matrices <- function(B, K, n.matrices = 10,
                                  alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
                                  alpha = 1.0, maxit = 200, verbose = F) {
  
  ### Inputs
  ## B: the matrix with the summary statistics
  ## K: the number of components
  ## n.matrices: the number of unmixing matrices to generate
  ## all other inputs are parameters for FastICA
  
  ### Output: a list of unmixing matrices
  
  m = nrow(B)
  t = ncol(B)
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = svd(B) 
  
  U = tsvd.list$u[,1:K]
  d = tsvd.list$d[1:K]
  V = tsvd.list$v[,1:K]
  
  G = cbind(t(U), t(V)) ## The components are in rows
  
  G.scaled = G*sqrt((m+t-1)/2) ## scaling to ensure that cov = I. It holds that m+t == nrow(G)
  
  unmixing.matrices = list()
  
  for (i in 1:n.matrices) {
    
    w.init = matrix(rnorm(K^2), K, K) ## using different initialization for each iteration
    
    if (alg.typ == "deflation") 
      unmixing.matrix = ica.R.def(X = G.scaled,
                                  n.comp = K,
                                  tol = tol, fun = fun,
                                  alpha = alpha, maxit = maxit, 
                                  verbose = verbose, w.init = w.init)
    
    else if (alg.typ == "parallel")
      unmixing.matrix = ica.R.par(X = G.scaled,
                                  n.comp = K,
                                  tol = tol, fun = fun,
                                  alpha = alpha, maxit = maxit,
                                  verbose = verbose, w.init = w.init)
    
    
    unmixing.matrices[[i]] = t(unmixing.matrix) ## transposed unmixing matrix, the components are in the columns
  }
  
  return(unmixing.matrices)
}



## Function to get a number of latent components by running ICA many times, as proposed by
## Lazarev et al. in the paper "GUIDE deconstructs genetic architectures using association studies",
## where GUIDE was introduced.

## Using different initializations with the same input, a pre-specified number of unmixing matrices
## is generated. Then, all possible pairs are compared to in terms of how many matching columns they have
## and the median is reported as the final estimate of the number of latent components.


get.nlatents <- function(B, starting.K, validation.reps = 10, cor.thres = 0.95,
                         alg.typ = "parallel", return.estimates = F) {
  
  ### Inputs
  ## B: the matrix with the summary statistics
  ## starting.K: the starting number of components
  ## validation.reps: the number of unmixing matrices to generate and compare with one another
  ## cor.thres: the correlation threshold (in absolute value) above which two columns of unmixing matrices are considered to match
  ## alg.typ: method for running FastICA, either "parallel" (components estimated at once) or "deflation" (components estimated sequentially)
  ## return.estimates: whether to return the estimated number of components from each pair of unmixing matrices comparison
  
  ### Output: a list containing the either the estimated number of latent components and the estimate's standard deviation,
  ###         or a vector of all individual estimated from pairwise comparisons

  unmixing.matrices = get.unmixing.matrices(B, K = starting.K, n.matrices = validation.reps, alg.typ = alg.typ) ## use the default parameters for FastICA
  
  combs = combn(1:validation.reps, 2)
  
  n.selected.components = c()
  
  for (i in 1:ncol(combs)) {
    
    idx1 = combs[1, i]
    idx2 = combs[2, i]
    
    A.1 = unmixing.matrices[[idx1]]
    A.2 = unmixing.matrices[[idx2]]
    
    components = 0
    
    for (ii in 1:starting.K) {
      for (jj in 1:starting.K) {
        
        if (abs(cor(A.1[,ii], A.2[,jj])) > cor.thres) { 
          
          components = components + 1
          break ## if matching columns are found stop searching and continue to the next column
        }
      }
    }
    n.selected.components = c(n.selected.components, components)
  }
  
  if (return.estimates) return(n.selected.components) ## the estimate from individual pairwise comparisons can be returned
  
  est.comp = median(n.selected.components)
  est.comp.sd = sd(n.selected.components)
  
  return(list(K = est.comp, sdev = est.comp.sd))
}


## The following two functions contain the implementation of the novel extension of GUIDE with ICASSO, 
## where the unmixing matrix is estimated  by incorporating many ICA runs. Details can be found in the thesis pdf.


get_optimal_unmixing_matrix <- function(unmix.total, clusters, cors.unmix) {
  
  ### See below for the description of inputs
  ### Output: The final unmixing matrix (transposed), estimated with ICASSO.
  
  K = length(unique(clusters))
  
  final.unmixing.matrix = matrix(0, K, K)
  
  for (ii in 1:K) {
    cluster_points = which(clusters == ii)  
    sub_cors = cors.unmix[cluster_points, cluster_points]  
    if (length(cluster_points) == 1) medoid_index = cluster_points
    else medoid_index = cluster_points[which.max(rowSums(sub_cors))]  
    final.unmixing.matrix[,ii] = unmix.total[,medoid_index]
  }
  return(final.unmixing.matrix)
}



get_guide <- function(B, K=10, ica_runs = 1, alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
                      alpha = 1.0, maxit = 200, verbose = T) {
  
  ### Inputs
  ## B: the matrix with the summary statistics
  ## K: the number of components
  ## ica_runs: the number of ICA runs used in ICASSO
  ## all other inputs are parameters for FastICA
  
  ### Output: a list consisting of:
  ## W.xl: The variant-to-latent GUIDE weights
  ## W.lt: The latent-to-trait GUIDE weights
  ## A.T:  The final unmixing matrix (transposed), estimated with ICASSO.
  ## hc: the hierarchical clustering object obtained from clustering all estimated components across ICA runs,
  ##     based on their absolute Pearson correlations.
  ## clusters: a vector of cluster assignments for each estimated component (length = K * ica_runs),
  ##           where each cluster corresponds to one of the K independent components.
  ## cors.unmix: the absolute correlation matrix of all components from all ICA runs,
  ##             used to compute similarity between estimated components.
  
  
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = svd(B)
  
  U = tsvd.list$u[,1:K]
  d = tsvd.list$d[1:K]
  V = tsvd.list$v[,1:K]
  
  
  if (ica_runs == 1) {
    
    A.T = get.unmixing.matrices(B, K = K, n.matrices = 1, verbose = F)[[1]]
    
    hc = NULL
    clusters = NULL
    cors.unmix = NULL
  }
  
  else {
    
    unmixing.matrices = get.unmixing.matrices(B, K = K, n.matrices = ica_runs, verbose = F) ## each component corresponds to a column in the unmixing matrix
    
    unmix.total = do.call(cbind, unmixing.matrices)
    
    cors.unmix = abs(cor(unmix.total))
    
    dist.unmix = 1 - cors.unmix ## distance measure between two components, defined as 1 minus their absolute Pearson correlation
    
    hc = hclust(as.dist(dist.unmix), method = "average")
    
    clusters = cutree(hc, k = K)
    
    A.T = get_optimal_unmixing_matrix(unmix.total, clusters, cors.unmix)
  }
  
  W.xl = U %*% A.T
  W.lt = V %*% A.T
  
  return(list(W.xl = W.xl,
              W.lt = W.lt,
              A=A.T,
              hc = hc,
              clusters = clusters,
              cors.unmix = cors.unmix))
} 
  
  
## Function to estimate the Cluster Quality Index (CQI) for each component.
## The CQI is defined as (within-cluster similarity) minus (between-cluster similarity).
## Higher values indicate more stable and well-separated components.


get_cqi_values <- function(cors.unmix, clusters) {
  
  ### See above for the inputs
  
  ### Output: a vector of CQI values
  
  cqi_values = c()
  K = length(unique(clusters))
  
  for (ii in 1:K) {
    
    in_cluster = which(clusters == ii)
    out_cluster = which(clusters != ii)
    
    within_sim = sum(cors.unmix[in_cluster, in_cluster])/(length(in_cluster)^2)
    between_sim = sum(cors.unmix[in_cluster, out_cluster])/(length(in_cluster)*length(out_cluster))
    
    cqi_values = c(cqi_values, within_sim - between_sim)
  }
  
  return(cqi_values)
}


## Function to compute the R index, a scalar value measuring overall clustering quality for different partitions,
## where lower values suggest more distinct and compact clusters. See  Levine and Domany (2001) for further details.


get_r_index <- function(cors.unmix, hc, K) {
  
  ### Inputs:
  ## cors.unmix, hc: as defined above
  ## K: the number of components
  
  ### Output: the R index
  
  
  dist.unmix = 1 - cors.unmix
  
  clusters = cutree(hc, K)
  
  I.R = 0
  
  for (ii in 1:K) {
    
    in_cluster = which(clusters == ii)
    within_dif = sum(dist.unmix[in_cluster, in_cluster])/(length(in_cluster)^2)
    
    out_indices = setdiff(1:K, ii)
    out_difs = c()
    
    for (m in out_indices) {
      
      out_cluster = which(clusters == m)
      out_difs = c(out_difs, sum(dist.unmix[in_cluster, out_cluster])/(length(in_cluster)*length(out_cluster)))
    }
    
    I.R = I.R + within_dif/min(out_difs)
  }
  return(I.R/K)
}




