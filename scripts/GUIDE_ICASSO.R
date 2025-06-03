
library(fastICA) ## package for implementation of fastICA


## Functions for the decomposition of a GWAS summary statistics matrix using the method
## Genetic Unmixing by Independent Decomposition (GUIDE).

## Function to get unmixing matrices from multiple ICA runs

get.unmixing.matrices <- function(B, K, n.matrices = 10,
                                  alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
                                  alpha = 1.0, maxit = 200, verbose = T) {
  
  m = nrow(B)
  t = ncol(B)
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = svd(B) #get_tsvd(B, K = starting.K)
  
  U = tsvd.list$u[,1:K]
  d = tsvd.list$d[1:K]
  V = tsvd.list$v[,1:K]
  
  G = cbind(t(U), t(V)) ## The components are in rows
  
  G.scaled = G*sqrt((m+t-1)/2) ## scaling to ensure that cov = I. It holds that m+t == nrow(G)
  
  unmixing.matrices = list()
  
  start_time = Sys.time()
  
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
  
  end_time = Sys.time()
  elapsed_time = as.numeric(difftime(end_time, start_time, units = "mins"))
  #print(paste(n.matrices, "unmixing matrices generated in", round(elapsed_time,2), "minutes"))
  
  return(unmixing.matrices)
}



## Function to get a number of latent components by running ICA many times, as proposed by
## Lazarev et al. in the paper "GUIDE deconstructs genetic architectures using association studies",
## where GUIDE was introduced.


get.nlatents <- function(B, starting.K, validation.reps = 10, cor.thres = 0.95,
                         alg.typ = "parallel", return.estimates = F) {
  
  ### Inputs
  ## B: the matrix with the summary statistics
  ## starting.K: the starting number of components
  ## validation.reps: the number of unmixing matrices to generate
  ## cor.thres: the threshold for which two columns of unmixing matrices are considered to match
  ## alg.typ, tol, fun, alpha, maxit, verbose: parameters used for the FastICA algorithm
  
  ### Output: a list containing the estimated number of components and the estimate's standard deviation

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


## Function for the novel extension of GUIDE with ICASSO...

guide_icasso <- function(B, K, ica_runs = 20) {
  
  unmixing.matrices = get.unmixing.matrices(B, K = K, n.matrices = ica_runs, verbose = F) ## each component corresponds to a column in the unmixing matrix
  
  unmix.total = do.call(cbind, unmixing.matrices)
  
  cors.unmix = abs(cor(unmix.total))
  
  dist.unmix = 1 - cors.unmix ## distance measure between two components, defined as 1 minus their absolute Pearson correlation
  
  hc = hclust(as.dist(dist.unmix), method = "average")
  
  clusters = cutree(hc, k = K)
  
  optimal.unmixing.matrix = matrix(0, K, K)
  
  for (ii in 1:K) {
    cluster_points = which(clusters == ii)  
    sub_cors = cors.unmix[cluster_points, cluster_points]  
    if (length(cluster_points) == 1) medoid_index = cluster_points
    else medoid_index = cluster_points[which.max(rowSums(sub_cors))]  
    optimal.unmixing.matrix[,ii] = unmix.total[,medoid_index]
  }
  
  return(list(hc = hc,
              clusters = clusters,
              cors.unmix = cors.unmix,
              optimal.unmixing.matrix = optimal.unmixing.matrix))
}


get_optimal_unmixing_matrix <- function(unmix.total, clusters, cors.unmix) {
  
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
  
  



get_cqi_values <- function(cors.unmix, clusters) {
  
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



##### to be removed

get_ica_clustering <- function(W, K, reps = 20) {

unmixing.matrices = get.unmixing.matrices(W, K = K, n.matrices = reps, verbose = F) ## each component corresponds to a column in the unmixing matrix
unmix.total = do.call(cbind, unmixing.matrices)

cors.unmix = abs(cor(unmix.total))
dist.unmix = 1 - cors.unmix

hc = hclust(as.dist(dist.unmix), method = "average") ## hierarchical clustering with average linkage

clusters = cutree(hc, k = K) ## cutting the tree at K components

return(list(hc = hc, clusters = clusters, cors.unmix = cors.unmix, unmix.total = unmix.total))
}


