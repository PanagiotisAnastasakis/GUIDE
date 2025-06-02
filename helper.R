
#library(irlba)
library(fastICA) ## package for implementation of fastICA
library(Homo.sapiens)
library(dplyr)



######################################################
######################################################
######################################################
######                                          ######
######      Part 1: Decomposition Methods       ######
######                                          ######
######################################################
######################################################
######################################################






######################################################
#                                                    #
#                  Method 1 -> DeGAs                 #
#                                                    #
######################################################


## Approach 1) Proportion of variance explained: we require that the first K components explain p% of the variance

## Approach 2) Above average eigenvalues. The "average" is computed as the intertia (sum of eigenvalues) divided by the rank of the matrix


get_tsvd <- function(B, K=0, method="average_eig", p=0.7) {
  
  ### Inputs
  ## B: the (standardized) summary statistics matrix 
  ## K: the number of components to keep. K=0 implies that an estimation method should be used
  ## method: the method to use for choosing the number of components. Options are "average_eig",
  ##         where components with above average eigenvalues are kept and "var_explained", where
  ##         we increase the number of components until the desired percentage of variance is explained.
  ## p: the percentage of variance to be explained when choosing components. Used only for method="var_explained"
  
  ### Output: a list with the TSVD of the summary statistics matrix
  
  stopifnot("Enter a valid value for K" = K>=0 & K<=min(nrow(B), ncol(B)))
  stopifnot("Enter a valid percentage of variance to be explained" = p>=0 & p<=1)
  
  matrix.svd = svd(B)
  
  U = matrix.svd$u
  d = matrix.svd$d
  V = matrix.svd$v
  
  if (K==0) {
    
    if (method == "average_eig") {
      
      I = sum(d^2) ## inertia equals sum of squared singular values
      rank.B = sum(d>0) ## rank of a matrix is the number of non-zero singular values
      
      K = which.min(d^2 > I/rank.B) - 1 ## the singular values are ordered
    }
    
    if (method == "var_explained") {
      
      vars = cumsum(d^2 / sum(d^2))
      K = which.max(vars > p)
    }
  }
  
  return(list(U=U[,1:K], d=d[1:K], V=V[,1:K]))
}




######################################################
#                                                    #
#                  Method 2 -> GUIDE                 #
#                                                    #
######################################################


get.nlatents <- function(B, starting.K, validation.reps = 10, cor.thres = 0.95,
                         alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
                         alpha = 1.0, maxit = 200, verbose = T,
                         return.estimates = F) {
  
  ### Inputs
  ## B: the matrix with the summary statistics
  ## starting.K: the starting number of components
  ## validation.reps: the number of unmixing matrices to generate
  ## cor.thres: the threshold for which two columns of unmixing matrices are considered to match
  ## alg.typ, tol, fun, alpha, maxit, verbose: parameters used for the FastICA algorithm
  
  ### Output: a list containing the estimated number of components and the estimate's standard deviation
  
  m = nrow(B)
  t = ncol(B)
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = get_tsvd(B, K = starting.K)
  
  U = tsvd.list$U
  d = tsvd.list$d
  V = tsvd.list$V
  
  G = cbind(t(U), t(V)) ## The components are in rows
  
  G.preprocessed = G*sqrt((m+t-1)/2) ## scaling to ensure that cov = I. m+t == nrow(G)
  
  unmixing.matrices = list()
  
  start_time = Sys.time()
  
  for (i in 1:validation.reps) {
    
    w.init = matrix(rnorm(starting.K^2), starting.K, starting.K) ## using different initialization for each iteration
    
    if (alg.typ == "deflation") {
      a = ica.R.def(X = G.preprocessed,
                    n.comp = starting.K,
                    tol = tol, fun = fun,
                    alpha = alpha, maxit = maxit, 
                    verbose = verbose, w.init = w.init)
    }
    else if (alg.typ == "parallel") {
      a = ica.R.par(X = G.preprocessed,
                    n.comp = starting.K,
                    tol = tol, fun = fun,
                    alpha = alpha, maxit = maxit,
                    verbose = verbose, w.init = w.init)
    }
    
    unmixing.matrices[[i]] = t(a)
  }
  
  end_time = Sys.time()
  elapsed_time = as.numeric(difftime(end_time, start_time, units = "mins"))
  print(paste(validation.reps, "unmixing matrices generated in", round(elapsed_time,2), "minutes"))
  
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
  
  if (return.estimates) return(n.selected.components)
  
  est.comp = median(n.selected.components)
  est.comp.sd = sd(n.selected.components)
  
  return(list(K = est.comp, sdev = est.comp.sd))
}



get.unmixing.matrices <- function(B, K, validation.reps = 10,
         alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
         alpha = 1.0, maxit = 200, verbose = T) {
  
  m = nrow(B)
  t = ncol(B)
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = get_tsvd(B, K = K)
  
  U = tsvd.list$U
  d = tsvd.list$d
  V = tsvd.list$V
  
  G = cbind(t(U), t(V)) ## The components are in rows
  
  G.preprocessed = G*sqrt((m+t-1)/2) ## scaling to ensure that cov = I. It holds that m+t == nrow(G)
  
  unmixing.matrices = list()
  
  start_time = Sys.time()
  
  for (i in 1:validation.reps) {
    
    w.init = matrix(rnorm(K^2), K, K) ## using different initialization for each iteration
    
    if (alg.typ == "deflation") 
      unmixing.matrix = ica.R.def(X = G.preprocessed,
                                  n.comp = K,
                                  tol = tol, fun = fun,
                                  alpha = alpha, maxit = maxit, 
                                  verbose = verbose, w.init = w.init)
    
    else if (alg.typ == "parallel")
      unmixing.matrix = ica.R.par(X = G.preprocessed,
                                  n.comp = K,
                                  tol = tol, fun = fun,
                                  alpha = alpha, maxit = maxit,
                                  verbose = verbose, w.init = w.init)
    
    
    unmixing.matrices[[i]] = t(unmixing.matrix)
  }
  
  end_time = Sys.time()
  elapsed_time = as.numeric(difftime(end_time, start_time, units = "mins"))
  #print(paste(validation.reps, "unmixing matrices generated in", round(elapsed_time,2), "minutes"))
  
  return(unmixing.matrices)
}



get_ica_clustering <- function(W, K, reps = 20) {
  
  unmixing.matrices = get.unmixing.matrices(W, K = K, validation.reps = reps, verbose = F) ## each component corresponds to a column in the unmixing matrix
  
  unmix.total = do.call(cbind, unmixing.matrices)
  
  cors.unmix = abs(cor(unmix.total))
  
  dist.unmix = 1 - cors.unmix
  
  hc = hclust(as.dist(dist.unmix), method = "average")
  
  clusters = cutree(hc, k = K)
  
  return(list(hc = hc, clusters = clusters, cors.unmix = cors.unmix, unmix.total = unmix.total))
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





get_guide <- function(B, K=10, unmixing.matrix = NA, alg.typ = "parallel", tol = 1e-04, fun = "logcosh",
                      alpha = 1.0, maxit = 200, verbose = T, use.irlba = F) {
  
  m = nrow(B)
  t = ncol(B)
  
  B = scale(B, center = TRUE, scale = FALSE)
  B = t(scale(t(B), center = TRUE, scale = FALSE))
  
  tsvd.list = get_tsvd(B, K = K, use.irlba = use.irlba)
  
  U = tsvd.list$U
  d = tsvd.list$d
  V = tsvd.list$V
  
  G = cbind(t(U), t(V)) ## The components are in rows
  
  G.preprocessed = G*sqrt((m+t-1)/2)
  
  w.init = matrix(rnorm(K^2), K, K) 

  if (!is.matrix(unmixing.matrix)) {
    
    print("Unmixing matrix generated")
    
    if (alg.typ == "deflation") 
       unmixing.matrix = ica.R.def(X = G.preprocessed,
                                   n.comp = K,
                                   tol = tol, fun = fun,
                                   alpha = alpha, maxit = maxit, 
                                   verbose = verbose, w.init = w.init)

    else if (alg.typ == "parallel")
       unmixing.matrix = ica.R.par(X = G.preprocessed,
                                   n.comp = K,
                                   tol = tol, fun = fun,
                                   alpha = alpha, maxit = maxit,
                                   verbose = verbose, w.init = w.init)

    A.T = t(unmixing.matrix)
    #print(A)
  }
  
  else A.T = unmixing.matrix ## the unmixing matrix has been given, assumed to be already transposed
  
  
  W.xl = U %*% A.T
  W.lt = V %*% A.T
  
  #W.lt = t(A.T) %*% t(V) 
  #W.lt = t(W.lt),
  
  return(list(A=A.T,
              W.xl = W.xl,
              W.lt = W.lt,
              U=U,
              d=d,
              V=V))
}























######################################################
######################################################
######################################################
######                                          ######
######    Part 2: Simulation Helper Functions   ######
######                                          ######
######################################################
######################################################
######################################################



## Function to find the correct permutation and sign for a set of simulated/estimated weights.

get_matched_weights <- function(weights.true, weights.estimated, n.variants,
                                return.only.indices = FALSE) {
  
  ### Inputs
  ## weights.true: the simulated (true) weights to and from the latent factors, concatenated. It has dimension (n+p)xL,
  ##               where n is the number of variants, p is the number of traits and L is the true number of latents.
  ## weights.estimated: the estimated weights to and from the latent factors concatenated, with dimension (n+p)xK,
  ##               where K is the number of components we estimated and can differ from L.
  
  ### Output: a matrix with two columns and min(L,K), the first containing the ordering of the estimated weights to
  ##          match the true ones and the second containing either 1 or -1, depending on whether the matched weights are
  ##          positively or negatively correlated.
  
  L = ncol(weights.true) 
  K = ncol(weights.estimated) 
  
  ## Capturing the correlation of each true set of weights and the estimated ones:
  
  cors.matrix = matrix(NA, L, K)
  for (ii in 1:L) {
    cors = apply(weights.estimated, 2, function(column) cor(weights.true[,ii], column))
    cors.matrix[ii,] = cors
  }
  
  
  matches = matrix(NA, min(L,K), 3) ## column 1 contains the indices of the true weights, column 2 is for the estimated weights indices
  ## and column 3 has -1 if their correlation is negative, otherwise has 1.
  
  for (ii in 1:min(L,K)) {
    
    max_idx = which(abs(cors.matrix) == max(abs(cors.matrix)), arr.ind = TRUE)
    i = max_idx[1] ## the true component
    j = max_idx[2] ## the estimated component
    
    if (cors.matrix[max_idx] < 0) matches[ii,3] = -1
    else matches[ii,3] = 1
    
    # Set all values in row i and column j to 0 to not consider them again
    cors.matrix[i, ] = 0
    cors.matrix[, j] = 0
    
    matches[ii, c(1,2)] = c(i,j)
    
  }
  
  if (return.only.indices) return(matches)
  
  if (K < L) {
    sorted.matches = matches[order(matches[, 2]), ]

    weights.matched = weights.true[,sorted.matches[,1]] %*% diag(sorted.matches[,3])
    
    weights.matched.xl = head(weights.matched, n.variants) ## we match the true ones !!
    weights.matched.lt = tail(weights.matched, nrow(weights.matched) - n.variants)
    
  }
  
  else {
    sorted.matches = matches[order(matches[, 1]),]
    
    weights.matched = weights.estimated[,sorted.matches[,2]] %*% diag(sorted.matches[,3]) ## we re-order the estimated components and multiply them with either -1 or 1 to correct for sign change
    
    weights.matched.xl = head(weights.matched, n.variants) ## we keep the first n rows, as these correspond to variants (W.xl.estimated) and the rest (below) are the traits (W.lt.estimated)
    weights.matched.lt = tail(weights.matched, nrow(weights.matched) - n.variants)
    
  }
  
  return(list(weights.matched.xl = weights.matched.xl, weights.matched.lt = weights.matched.lt))
}





## Function to plot either the variants-to-latents or the latents-to-traits weights, both assumed to be matched.



plot_simulation <- function(weights.true, weights.estimated, xlab, ylab, main) {
  
  ### Inputs
  ## weights.true: the simulated (true) weights to or from the latent factors
  ## weights.estimated: the estimated weights to or from the latent factors
  ## xlab, ylab, main: used for plotting
  
  ### Output: a scatter plot of simulated against true weights

  L <- ncol(weights.true)  
  df <- data.frame(
    True = as.vector(weights.true),
    Estimated = as.vector(weights.estimated),
    Component = rep(1:L, each = nrow(weights.true))
  )
  
  colors <- rainbow(L)
  
  ggplot(df, aes(x = Estimated, y = True, color = factor(Component))) +
    geom_point(alpha = 0.7) + 
    scale_color_manual(values = colors) + 
    labs(title = main, x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 25, hjust = 0.5, vjust = 3),
          axis.title.x = element_text(size = 22, vjust = -3),
          axis.title.y = element_text(size = 22, vjust = 3),
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 16),
          plot.margin = margin(20, 15, 25, 15))
}





### Function to generate latent associations with a given level of dependency

generate_latent_effect_matrix <- function(n, n.latents, lower, upper, dependence.prob = 0.1, noise.sd = 0.02) {
  
  ### Inputs
  ## n: the number of vertices (either variants or traits)
  ## n.latents: the number of latents of the simulated structure
  ## lower, upper: the lower and upper bounds for the Uniform distribution from which the weights are sampled
  ## dependence.prob: for this simulation strategy, we initially assume that the vertices are clustered in `n.latents` different clusters of same size `n`/`n.latents`,
  ##                  and that the effects to/from the i'th latent are sampled from the vector `effs` for the i'th cluster and are 0 for all other clusters.
  ##                  Because this scenario is idealized, we add the option of introducing dependencies, but assigning a probability
  ##                  `dependence.prob` for each vertex outside the cluster to have an effect, randomly sampled from the vector `effs`.
  ## noise.sd: the standard deviation for sampling a matrix from the gaussian distribution with mean 0 to add to the generated latent effect matrix
  
  ### Output: the latent effect matrix (either Wxl or Wlt)
  
  W = matrix(0, nrow = n, ncol = n.latents)
  
  for (i in seq_len(n.latents)) {
    
    row_start = (i - 1) * (n / n.latents) + 1
    row_end = i * (n / n.latents)
    
    #block = sample(effs, (n / n.latents), replace=T)
    block = runif(n / n.latents, min = lower, max = upper)
    
    W[row_start:row_end, i] = block
    
    if (dependence.prob > 0) {
      remaining.effs = replicate(n - n/n.latents, {
        if (runif(1) < dependence.prob) {
          #sample(effs, 1)
          runif(1, min = lower, max = upper)
        } else {
          0
        }
      })
      
      W[-(row_start:row_end), i] = remaining.effs
    }
  }
  
  W = W + matrix(rnorm(n*n.latents, 0, noise.sd), nrow = n, ncol = n.latents) ## add noise
  
  return(W)
}



## Function for simulations under a polygenic additive model


polygenic_sim <- function(n, p, k, sampler = "gaussian", sparsity = "none", h.sq = 0.1, Nsamples = 1e5, return.afs = F) {
  
  ## Nsamples controls the power and thus the standard error. Higher sample size means smaller SE
  
  afs = runif(n, 0.01, 0.5)
  
  xl.variance = h.sq/(n*k*2*(afs)*(1-afs))
  lt.variance = 1
  
  sparsity = switch(sparsity,
                    "none" = rep(0,k),
                    "low" = runif(k, 0, 0.2),
                    "high" = runif(k, 0.2, 0.7))
  
  if (sampler == "gaussian") {
    
    W.xl = do.call("cbind", lapply(sparsity, function(s) rnorm(n, sd=sqrt(xl.variance)) * rbinom(n,size = 1,prob = 1-s)))
    
    while (any(apply(W.xl,2,function(x) sum(x==0)) == nrow(W.xl))) { ## require none of the columns are 0
      W.xl = do.call("cbind", lapply(sparsity, function(s) rnorm(n, sd=sqrt(xl.variance)) * rbinom(n,size = 1,prob = 1-s)))
    }
    
    W.lt = do.call("cbind", lapply(sparsity, function(s) rnorm(p, sd=sqrt(lt.variance)) * rbinom(p,size = 1,prob = 1-s)))
    
    while (any(apply(W.lt,2,function(x) sum(x==0)) == nrow(W.lt))) { ## require none of the columns are 0
      W.lt = do.call("cbind", lapply(sparsity, function(s) rnorm(p, sd=sqrt(lt.variance)) * rbinom(p,size = 1,prob = 1-s)))
    }
  }
  
  else if (sampler == "laplace") {
    
    W.xl = do.call("cbind", lapply(sparsity, function(s) rLaplace(n, b=sqrt(xl.variance/2)) * rbinom(n,size = 1,prob = 1-s)))
    
    while (any(apply(W.xl,2,function(x) sum(x==0)) == nrow(W.xl))) { ## require none of the columns are 0
      W.xl = do.call("cbind", lapply(sparsity, function(s) rLaplace(n, b=sqrt(xl.variance/2)) * rbinom(n,size = 1,prob = 1-s)))
    }
    
    W.lt = do.call("cbind", lapply(sparsity, function(s) rLaplace(p, b=sqrt(lt.variance/2)) * rbinom(p,size = 1,prob = 1-s)))
    
    while (any(apply(W.lt,2,function(x) sum(x==0)) == nrow(W.lt))) { ## require none of the columns are 0
      W.lt = do.call("cbind", lapply(sparsity, function(s) rLaplace(p, b=sqrt(lt.variance/2)) * rbinom(p,size = 1,prob = 1-s)))
    }
  }
  
  B = W.xl %*% t(W.lt)
  
  SEs = sqrt(1/(Nsamples*2*afs*(1-afs)))
  
  E = do.call("cbind", lapply(1:p, function(j) rnorm(n, sd = SEs))) ##noise (assuming independent SNPs)
  
  B = B + E ## adding noise
  
  W = B / replicate(p, SEs) ## dividing by the standard errors to get the z-scores
  
  if (return.afs) return(list(W = W, W.xl = W.xl, W.lt = W.lt, afs = afs))
  
  return(list(W = W, W.xl = W.xl, W.lt = W.lt))
}











######################################################
######################################################
######################################################
######                                          ######
######  Part 3: Data Analysis Helper Functions  ######
######                                          ######
######################################################
######################################################
######################################################



### Important note here::: with GUIDE/DeGAs, we are able to quantify the contribution of different SNPs/phenotypes in terms of magnitude,
### but we are unable to know whether their effect is positive or negative (as highlighted in FactorGo) due to possible sign changes.


# contr_genes <- function(W.xl, genes, K) {
#   
#   ### Inputs
#   ## U: the variant matrix resulting from TSVD
#   ## genes: the indices used to compute the contribution scores. If it is a single value, it corresponds to a variant.
#   ##        If it is a vector of variants, it corresponds to a gene. Otherwise, it is a list of genes.
#   ## K: the component for which we compute the contribution.
#   
#   ### Output: the contribution score for each gene
#   
#   if (is.list(genes)) {
#     
#     common_variants = Reduce(intersect, genes)
#     
#     return(sum(W.xl[common_variants,K]^2))
#   }
#   
#   return(W.xl[genes,K]^2)
# }
# 
# 
# contr_phenotype <- function(W.lt, phenotype, K) {
#   
#   ### Inputs
#   ## V: the phenotype matrix resulting from TSVD
#   ## phenotype: the index for the phenotype
#   ## K: the component for which we compute the contribution
#   
#   ### Output: the contribution score
#   
#   return(W.lt[phenotype,K]^2)
# }




# sq_cosine_scores <- function(tsvd.list, K, index, variant=T) {
#   
#   ### Inputs
#   ## tsvd.list: the TSVD of the summary statistics matrix (a list with U,D,V)
#   ## K: the component
#   ## index: the index of the phenotype or variant for which we compute the contribution
#   ## variant: True if we compute the contribution of a variant, false if we compute the contribution of a phenotype
#   
#   ### Output: the contribution score
#   
#   U = tsvd.list$U
#   d = tsvd.list$d
#   V = tsvd.list$V
#   
#   if (variant) {
#     
#     U = tsvd.list$U
#     d = tsvd.list$d
#     
#     factor.score = U%*%diag(d)
#   }
#   
#   else {
#     
#     V = tsvd.list$V
#     d = tsvd.list$d
#     
#     factor.score = V%*%diag(d)
#   }
#   
#   return(factor.score[index, K]^2 / sum(factor.score[index,]^2))
# }



#contr_genes(tsvd.list$U, genes = 20, K = 4) ## how much does variant 20 contribute to component 4?
#contr_genes(tsvd.list$U, genes = c(20,30,40), K = 4) ## how much does gene (20,30,40) contribute to component 4?
#contr_genes(tsvd.list$U, genes = list(c(5, 10, 20, 30, 40), c(10, 20, 25, 30, 35, 40), c(2, 3, 4, 20, 22, 30, 39, 40, 50)), K = 4) ## how much do genes ... contribute to component 4?

#contr_phenotype(tsvd.list$V, phenotype=20, K= 2) ## how much does phenotype 20 contribute to component 2?

#sq_cosine_scores(tsvd.list, K=1, 13, variant=F)  ## How much does component 1 contribute to phenotype 13?
#sq_cosine_scores(tsvd.list, K=3, 3, variant=T)  ## How much does component 3 contribute to variant 3?


####### Source: https://rdrr.io/github/pws3141/clusterICA/src/R/entropy.R

#'Calculates entropy using m-spacing.
#'
#'
#' @param x the data, either a vector or matrix.
#'              If x is a matrix, entropy is estimated for each row separately.
#' @param m (optional) the m-spacing. Defaults to m <- sqrt(n)
#'              if missing, where n is length(x) if x is a vector,
#'                  or ncol(x) if a matrix
#'
#' @return Vector of real numbers corresponding to the approximate
#'              entropy for each row of input x.

# mSpacingEntropy <- function(x, m) {
#   
#   if(is.vector(x)) x <- matrix(x, nrow=1)
#   if(ncol(x) == 1) stop("require p > 1")
#   xt <- apply(x, 1, function(x) sort(x)) # nb: change to xt
#   n <- nrow(xt)
#   if(missing(m)) m <- floor(sqrt(n))
#   
#   d <- xt[(m+1):n,, drop=FALSE] - xt[1:(n-m),, drop=FALSE]
#   apply(d, 2, function(dd) {
#     (1/n) * sum(log((n / m) * dd))
#   }) - digamma(m) + log(m)
# }


## Note that we sum the entropies for every x->l or l->t component. Is that correct??







# Genetics helper

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Homo.sapiens")
# BiocManager::install("GenomicRanges")


# Code from bNMF for T2D pipeline file format_bNMF_results.Rmd
# snps is a vector of variants in the form "chr:pos" or "chr:pos:ref:alt"
# Returns a data frame with variant name, gene name and ENTREZID. 
# If a gene name was not found, defaults to variant name

query_locus_names <- function(snps){
  
  geneRanges <- function(db, column="ENTREZID"){
    g <- genes(db, columns=column)
    col <- mcols(g)[[column]]
    genes <- granges(g)[rep(seq_along(g), elementNROWS(col))]
    mcols(genes)[[column]] <- as.character(unlist(col))
    genes
  }
  
  gns = geneRanges(Homo.sapiens, column="SYMBOL")
  # need chr | start | end, with gene names for row names
  df_gns <- data.frame(gns)
  gtf.gene <- df_gns %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    dplyr::select(chr, start, end, SYMBOL)
  
  # need chr | start | end, with gene names for row names
  df_entrez <- geneRanges(Homo.sapiens, column="ENTREZID") %>%
    data.frame() %>%
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>% 
    dplyr::select(ChrPos, ENTREZID)
  
  #' Convert from string to range
  #' 
  #' @param pos A vector of strings ex. chr1 2938302 2938329
  #' @param delim Delimiter for string splitting
  #' @param region Boolean of whether region or just one position
  #'
  #' @returns Dataframe of ranges
  #' 
  string2range <- function(pos, delim=' ', region=TRUE) {
    posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
    posp[,1] <- posp[,1]
    posp[,2] <- as.numeric(as.character(posp[,2]))
    if(region) {
      posp[,3] <- as.numeric(as.character(posp[,3]))
    } else {
      posp[,3] <- posp[,2]
    }
    return(posp)
  }
  
  #' Convert from ranges to GRanges
  #' 
  #' @param df Dataframe with columns as sequence name, start, and end
  #' 
  #' @returns GRanges version 
  #' 
  range2GRanges <- function(df) {
    require(GenomicRanges)
    require(IRanges)
    gr <- GenomicRanges::GRanges(
      seqnames = df[,1],
      ranges=IRanges(start = df[,2], end = df[,3])
    )
    return(gr)
  }
  
  # convert SNPs to GRanges
  snps.ranges <- string2range(snps, delim=":", region=FALSE)
  snps.granges <- range2GRanges(snps.ranges)
  names(snps.granges) <- snps
  
  # convert genes to GRanges
  gtf.granges <- range2GRanges(gtf.gene)
  names(gtf.granges) <-  gtf.gene$SYMBOL  #gene.names
  
  hits <- GenomicRanges::nearest(snps.granges, gtf.granges)
  # make vector of SNPs to gene
  
  df_gns <- df_gns %>% 
    mutate(chr = gsub("chr","",seqnames)) %>%
    mutate(ChrPos = paste(chr,start,sep=":")) %>%
    dplyr::select(chr, start, end, ChrPos, SYMBOL) %>%
    merge(df_entrez, by="ChrPos")
  
  df_hits <- data.frame(gene=names(gtf.granges)[hits]) %>%
    mutate(variant = names(snps.granges)) %>%
    merge(df_gns[,c("SYMBOL","ENTREZID")], by.x="gene", by.y="SYMBOL") %>%
    filter(!duplicated(variant))
  
  # Make the duplicated genes have unique names.
  #If no gene is found for a variant, then the name of the variant is returned
  if (nrow(df_hits) > 0) {
    dup_genes <- df_hits$gene[!is.na(df_hits$gene)]
    duplicates <- data.frame(duplicated = table(dup_genes) > 1)
    duplicates$gene <- rownames(duplicates)
    
    for (r in 1:nrow(duplicates)) {
      row <- duplicates[r,]
      if (!is.na(row$duplicated) && row$duplicated) {
        idx <- which(df_hits$gene == row$gene)
        df_hits$gene[idx] <- paste0(row$gene, "_", seq_along(idx))
      }
    }
  }
  
  # else{ ## the default
  #   # Make the duplicated genes have unique names
  #   duplicates <- data.frame(duplicated = table(df_hits$gene) > 1)
  #   duplicates$gene <- rownames(duplicates)
  #   for(r in 1:nrow(duplicates)){
  #     row <- duplicates[r,]
  #     if(row$duplicated){
  #       df_hits$gene[df_hits$gene == row$gene] <- paste0(row$gene, "_", 1:sum(df_hits$gene == row$gene))
  #     }
  #   }
  # }
  # 
  
  res <- data.frame(variant = snps) %>%
    left_join(df_hits, by="variant") %>%
    mutate(gene = dplyr::coalesce(gene, variant)) # coalese fills in empty genes w/ variant
  
  return(res)
}






########## Visualization



cluster_circle_plot <- function(cluster_weights, total_names, title, rotate = 0) {
  
  cluster_data = data.frame(
    names = total_names,
    group = c(ifelse(cluster_weights[1:650] > 0, "Vpos", "Vneg"), ifelse(cluster_weights[651:760] > 0, "Tpos", "Tneg")),
    weights = cluster_weights*35
  )
  
  cluster_data = cluster_data[cluster_data$weights != 0, ]
  cluster_data = cluster_data %>% arrange(group, abs(weights))
  cluster_data = cluster_data %>% arrange(desc(abs(weights))) %>% head(30)
  
  empty_bar <- 4
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(factor(cluster_data$group)), ncol(cluster_data)))
  colnames(to_add) <- colnames(cluster_data)
  to_add$group <- rep(levels(factor(cluster_data$group)), each=empty_bar)
  cluster_data <- rbind(cluster_data, to_add)
  cluster_data <- cluster_data %>% arrange(group)
  cluster_data$id <- seq(1, nrow(cluster_data))
  
  rot = 0.5
  p.rot = 14.8
  if (nrow(cluster_data) < 30) {
    rot = 0.2
    p.rot = 5.27
  }
  if (length(unique(cluster_data$group)) < 4 & nrow(cluster_data) == 30) {
    rot = 0.2
    p.rot = 13.6
  }
  
  label_data <- cluster_data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - rot) / number_of_bar - rotate
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  group_colors <- c("Vpos" = "#228B22", "Vneg" = "#9370DB",
                    "Tpos" = "#F08080", "Tneg" = "#87CEEB")
  
  pl = ggplot(cluster_data, aes(x = as.factor(id), y = weights, fill = group)) +
    geom_bar(stat = "identity", alpha = 0.5) +
    ylim(-50, 50) +                  
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      #plot.margin = unit(rep(-1, 4), "cm"),
      plot.title = element_text(hjust = 0.5, size = 40, face = "bold", vjust=-11),
      plot.margin = margin(t = -70, r = -60, b = -20, l = -100, unit = "pt")   # <-- HERE: increase top (t)
    ) +
    coord_polar(start = rotate * pi/180) +
    scale_fill_manual(values = group_colors) +
    geom_text(data = label_data, 
              aes(x = id, 
                  y = ifelse(weights > 0, weights + 5, 5),
                  label = names, 
                  hjust = hjust), 
              color = "black", fontface = "bold", alpha = 0.6, 
              angle = label_data$angle, 
              inherit.aes = FALSE,
              size = 7) +
    annotate("path", x = seq(0, p.rot*pi, length.out = 200), y = rep(0, 200), 
             linewidth = 0.3, color = "black") +
    #annotate("text", x = 0, y = 50, label = title, size = 10, fontface = "bold")
    ggtitle(title)
  
  suppressWarnings(print(pl))
}















top_n_values_indices <- function(vec, n) {
  indices <- order(vec, decreasing = TRUE)[1:n]  # Get top n indices
  values <- vec[indices]  # Get corresponding values
  return(list(values = values, indices = indices))
}




stacked_barplot <- function(contr.scores, n, main = NULL, ylab = NULL, xlab = NULL) {
  
  k = ncol(contr.scores)
  
  data.matrix = matrix(0, nrow = n+1, ncol = k)

  for (ii in 1:k) {
         
         contr.lat = top_n_values_indices(contr.scores[,ii], n=n)
         data.matrix[,ii] = c(contr.lat$values, 1 - sum(contr.lat$values))
  }
           
  colors = c(colorRampPalette(c("darkblue", "blue", "cyan"))(n), "grey")
  
  default.par = par(no.readonly = TRUE)  
  par(mar = c(4, 4, 4, 2))  # Increase bottom margin to fit rotated labels
  
  bar_positions = barplot(data.matrix, beside = FALSE, col = colors, ylim = c(0, 1.05),
                          ylab = "", yaxt = "n", yaxs = "i", names.arg = rep("", k))
  
  axis(2, at = seq(0, 1, by = 0.2))
  mtext(ylab, side = 2, line = 2.5, adj = 0.475, cex = 1.4)
  title(main = main, cex.main = 2)
  
  labels <- paste0(xlab, 1:k)
  text(
    x = bar_positions,
    y = par("usr")[3] - 0.03,  # slightly below x-axis
    labels = labels,
    srt = 60,
    adj = 1,
    xpd = TRUE,
    cex = 0.9
  )
  
  par(default.par)
}




### Weighted Jaccard Index Computation

weighted_jaccard <- function(a, b) {
  
  jaccard.index = sum(pmin(a, b)) / sum(pmax(a, b))
  
  return(jaccard.index)
}

compute_jaccard_matrix <- function(A, B) {
  
  k.1 = ncol(A)
  k.2 = ncol(B)
  
  jaccard.matrix = matrix(0, nrow = k.1, ncol = k.2)
  
  for (ii in 1:k.1) {
    for (jj in 1:k.2) {
      jaccard.matrix[ii, jj] <- sum(pmin(A[,ii], B[,jj])) / sum(pmax(A[,ii], B[,jj]))
    }
  }
  
  return(jaccard.matrix)
}















####################################################################
####################################################################
####################################################################
####################################################################
####################################################################




#### REDUNDANT ####




######################################################
#                                                    #
#                  Method 3 -> bNMF                  #
#                                                    #
######################################################


## Function for NN-transformation


NN_transformation <- function(z_mat) {
  
  z_mat_pos <- z_mat
  z_mat_pos[z_mat_pos < 0] <- 0
  colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
  z_mat_neg <- -z_mat
  z_mat_neg[z_mat_neg < 0] <- 0
  colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
  final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  return(final_z_mat)
}



## bNMF code from: https://github.com/gwas-partitioning/bnmf-clustering


BayesNMF.L2EU <- function(V0, n.iter=10000, a0=10, tol=1e-7, K=15, K0=15, phi=1.0) {
  
  # Bayesian NMF with half-normal priors for W and H
  # V0: input z-score matrix (variants x traits)
  # n.iter: Number of iterations for parameter optimization
  # a0: Hyper-parameter for inverse gamma prior on ARD relevance weights
  # tol: Tolerance for convergence of fitting procedure
  # K: Number of clusters to be initialized (algorithm may drive some to zero)
  # K0: Used for setting b0 (lambda prior hyper-parameter) -- should be equal to K
  # phi: Scaling parameter
  
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[, active_nodes]
  V <- V0 - min(V0)
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  W <- matrix(runif(N * K) * Vmax, ncol=K)
  H <- matrix(runif(M * K) * Vmax, ncol=M)
  
  I <- array(1, dim=c(N, M))
  V.ap <- W %*% H + eps
  
  phi <- sd(V)^2 * phi
  C <- (N + M) / 2 + a0 + 1
  b0 <- 3.14 * (a0 - 1) * mean(V) / (2 * K0)
  lambda.bound <- b0 / C
  lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
  lambda.cut <- lambda.bound * 1.5
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while (del >= tol & iter < n.iter) {
    H <- H * (t(W) %*% V) / (t(W) %*% V.ap + phi * H * matrix(rep(1 / lambda, M), ncol=M) + eps)
    
    V.ap <- W %*% H + eps
    
    W <- W * (V %*% t(H)) / (V.ap %*% t(H) + phi * W * t(matrix(rep(1 / lambda, N), ncol=N)) + eps)
    
    V.ap <- W %*% H + eps
    
    lambda <- (0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / C
    
    del <- max(abs(lambda - n.lambda[[iter - 1]]) / n.lambda[[iter - 1]])
    
    like <- sum((V - V.ap)^2) / 2
    
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi * sum((0.5 * colSums(W^2) + 0.5 * rowSums(H^2) + b0) / lambda + C * log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V - V.ap)^2)
    
    if (iter %% 100 == 0) {
      cat(iter, n.evid[[iter]], n.like[[iter]], n.error[[iter]], del, 
          sum(colSums(W) != 0), sum(lambda >= lambda.cut), '\n')
    }
    iter <- iter + 1
  }
  return(list(
    W = W,  # Variant weight matrix (N x K)
    H = H,  # Trait weight matrix (K x M)
    n.like = n.like,  # List of reconstruction errors (sum of squared errors / 2) per iteration
    n.evid = n.evid,  # List of negative log-likelihoods per iteration
    n.lambda = n.lambda,  # List of lambda vectors (shared weights for each of K clusters, some ~0) per iteration
    n.error = n.error  # List of reconstruction errors (sum of squared errors) per iteration
  ))
}



run_bNMF <- function(z_mat, n_reps=10, random_seed=1, K=20, K0=10, tolerance=1e-7) {
  
  # Given an input matrix as created by prep_z_matrix(), run the bNMF procedure
  # a series of times to generate results and evaluate cluster stability
  
  set.seed(random_seed)
  
  bnmf_reps <- lapply(1:n_reps, function(r) {
    print(paste("ITERATION",r))
    res <- BayesNMF.L2EU(V0 = z_mat, a0 = log(sum(dim(z_mat))), K=K, K0=K0, tol=tolerance)
    names(res) <- c("W", "H", "n.like", "n.evid", "n.lambda", "n.error")
    res
  })
  return(bnmf_reps)
}


make_run_summary <- function(reps) {
  
  # Given a list of bNMF iteration outputs, summarize the K choices and associated likelihoods across runs
  
  run_summary <- map_dfr(1:length(reps), function(i) {
    res <- reps[[i]]
    final_lambdas <- res$n.lambda[[length(res$n.lambda)]]
    tibble(
      run=i,
      K=sum(final_lambdas > min(final_lambdas)),  # Assume that lambdas equal to the minimum lambda are ~ 0
      evid=res$n.evid[[length(res$n.evid)]]  # Evidence = -log_likelihood
    )
  }) %>%
    arrange(evid)
  
  unique.K <- table(run_summary$K)
  n.K <- length(unique.K)  # Number of distinct K
  MAP.K.run <- sapply(names(unique.K), function(k) {  # bNMF run index with the maximum posterior for given K
    tmp <- run_summary[run_summary$K == k, ]
    tmp$run[which.min(tmp$evid)]
  })
  
  list(run_tbl=run_summary, unique.K=unique.K, MAP.K.run=MAP.K.run)
}







# library(MASS)
# 
# 
# generate_data <- function(n.samples, n.features, var.explained, noise.level) {
#   
#   ### Inputs
#   ## n.samples: the number of samples to have in the simulated matrix
#   ## n.features: the number of samples to have in the simulated matrix
#   ## n.pcs: the number of significant principal components, i.e. those who explain more variance than the rest
#   ## noise.level: a scalar between 0 and 1 indicating the fraction of standard deviation of the data that is used as noise standard deviation
#   
#   ### Output: a simulated data matrix
#   
#   eigenvalues = diag(var.explained)
#   
#   rotation_matrix = qr.Q(qr(matrix(rnorm(n.features^2), n.features, n.features)))
#   
#   cov_matrix = rotation_matrix %*% eigenvalues %*% t(rotation_matrix)
#   
#   data = mvrnorm(n.samples, mu = rep(0, n.features), Sigma = cov_matrix)
#   
#   noise_sd = noise.level * sd(data)
#   
#   noise = matrix(rnorm(n.samples * n.features, mean = 0, sd = noise_sd), nrow = n.samples, ncol = n.features)
#   
#   data = data + noise
#   
#   return(data)
# }
# 
# 
# 
# generate_ve <- function(n.features, n.pcs, method) {
#   
#   eig = sort(rgamma(5*n.pcs, shape = 1), decreasing = T)[1:n.pcs]
#   
#   if (method == "uniform") {
#     
#     var.explained = c(eig, sort(runif(n.features - n.pcs, max = eig[n.pcs]/5), decreasing = T))
#   }
#   
#   if (method == "constant") {
#     
#     var.explained = c(eig, rep(eig[n.pcs]/30, n.features - n.pcs))
#   }
#   
#   if (method == "zero") {
#     
#     var.explained = c(eig, rep(0, n.features - n.pcs))
#   }
#   
#   var.explained = var.explained/sum(var.explained) ## scale so variance explained sums to 1
#   
#   return(var.explained)
# } 
# 
# 
# 
# 
# ## Function to generate block matrices
# 
# generate_block_matrix <- function(n, p, n.blocks, effs, noise.sd = 0.02) {
#   
#   ### Inputs
#   ## n: the number of rows
#   ## p: the number of columns
#   ## n.blocks: the number of blocks in the matrix
#   ## sds: the standard deviation of the samples in each block. If it is scalar, same sd for all blocks is assumed
#   ## noise.level: the percentage of the standard deviation of the initial matrix B that we use to generate noise, which we add to B
#   
#   ### Output: the block matrix B
#   
#   #if (length(sds) == 1) sds = rep(sds, n.blocks)
#   
#   W = matrix(0, nrow = n, ncol = p)
#   
#   
#   for (i in seq_len(n.blocks)) {
#     
#     row_start = (i - 1) * (n / n.blocks) + 1
#     row_end = i * (n / n.blocks)
#     col_start = (i - 1) * (p / n.blocks) + 1
#     col_end = i * (p / n.blocks)
#     
#     #block = matrix(rnorm((n / n.blocks) * (p / n.blocks), mus[i], 0.1),
#     #               n / n.blocks,
#     #               p / n.blocks)
#     
#     block = matrix(sample(effs, (n / n.blocks) * (p / n.blocks), replace=T),
#                    n / n.blocks,
#                    p / n.blocks)
#     
#     W[row_start:row_end, col_start:col_end] = block
#   }
#   
#   W = W + matrix(rnorm(n*p, 0, noise.sd), nrow = n, ncol = p) ## add noise
#   
#   return(W)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# permutation_test <- function(sample1, sample2, num.perm = 10000) {
#   
#   observed.stat = abs(mean(sample1) - mean(sample2))
#   
#   combined.samples = c(sample1, sample2)
#   n.1 = length(sample1) ## used in case that sample1, sample2 are not of equal length
#   
#   counts = 0
#   
#   for (i in 1:num.perm) {
#     
#     permuted.samples = sample(combined.samples)
#     
#     perm.sample1 = permuted.samples[1:n.1]
#     perm.sample2 = permuted.samples[(n.1 + 1):length(combined.samples)]
#     
#     permuted.stat = abs(mean(perm.sample1) - mean(perm.sample2))
#     
#     if (permuted.stat >= observed.stat) {
#       counts <- counts + 1
#     }
#   }
#   
#   p.value = counts / num.perm
#   
#   return(p.value)
# }
# 
# 
# 
# get_matched_weights_redundant <- function(weights.true, weights.estimated, n.variants) {
#   
#   ### Inputs
#   ## weights.true: the simulated (true) weights to and from the latent factors, concatenated. It has dimension (n+p)xL,
#   ##               where n is the number of variants, p is the number of traits and L is the true number of latents.
#   ## weights.estimated: the estimated weights to and from the latent factors concatenated, with dimension (n+p)xK,
#   ##               where K is the number of components we estimated and can differ from L.
#   
#   ### Output: a matrix with two columns and min(L,K), the first containing the ordering of the estimated weights to
#   ##          match the true ones and the second containing either 1 or -1, depending on whether the matched weights are
#   ##          positively or negatively correlated.
#   
#   L = ncol(weights.true) 
#   K = ncol(weights.estimated) 
#   
#   ## Capturing the correlation of each true set of weights and the estimated ones:
#   
#   cors.matrix = matrix(NA, L, K)
#   for (ii in 1:L) {
#     cors = apply(weights.estimated, 2, function(column) cor(weights.true[,ii], column))
#     cors.matrix[ii,] = cors}
#   
#   
#   matches = matrix(NA, min(L,K), 3) ## column 1 contains the indices of the true weights, column 2 is for the estimated weights indices
#   ## and column 3 has -1 if their correlation is negative, otherwise has 1.
#   
#   for (ii in 1:min(L,K)) {
#     
#     max_idx = which(abs(cors.matrix) == max(abs(cors.matrix)), arr.ind = TRUE)
#     i = max_idx[1] ## the true component
#     j = max_idx[2] ## the estimated component
#     
#     if (cors.matrix[max_idx] < 0) matches[ii,3] = -1
#     else matches[ii,3] = 1
#     
#     # Set all values in row i and column j to 0 to not consider them again
#     cors.matrix[i, ] = 0
#     cors.matrix[, j] = 0
#     
#     matches[ii, c(1,2)] = c(i,j)
#     
#   }
#   
#   #if (K < L) {
#   #  sorted.matches = matches[order(matches[, 2]), ]
#   #  matching.indices = sorted.matches[,c(1,3)]
#   #}
#   
#   sorted.matches = matches[order(matches[, 1]),] ## we sort the indices based on the ordering of the true components
#   matching.indices = sorted.matches[,c(2,3)] ## with find_corresponding_components, we were returning this, that's why we keep this line
#   
#   weights.matched = weights.estimated[,matching.indices[,1]] %*% diag(matching.indices[,2]) ## we re-order the estimated components and multiply them with either -1 or 1 to correct for sign change
#   
#   weights.matched.xl = head(weights.matched, n.variants) ## we keep the first n rows, as these correspond to variants (W.xl.estimated) and the rest (below) are the traits (W.lt.estimated)
#   weights.matched.lt = tail(weights.matched, nrow(weights.matched) - n.variants)
#   
#   return(list(weights.matched.xl = weights.matched.xl, weights.matched.lt = weights.matched.lt))
#   
# }



get_r_index <- function(cors.unmix, hc, L) {
  
  dist.unmix = 1 - cors.unmix
  
  clusters = cutree(hc, L)
  
  I.R = 0
  
  for (ii in 1:L) {
    
    in_cluster = which(clusters == ii)
    within_dif = sum(dist.unmix[in_cluster, in_cluster])/(length(in_cluster)^2)
    
    out_indices = setdiff(1:L, ii)
    out_difs = c()
    
    for (m in out_indices) {
      
      out_cluster = which(clusters == m)
      out_difs = c(out_difs, sum(dist.unmix[in_cluster, out_cluster])/(length(in_cluster)*length(out_cluster)))
    }
    
    I.R = I.R + within_dif/min(out_difs)
  }
  return(I.R/L)
}




library(grid)
library(gtable)


heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                           tree_row, treeheight_col, treeheight_row, filename, width, 
                           height, breaks, color, legend, annotation_row, annotation_col, 
                           annotation_colors, annotation_legend, annotation_names_row, 
                           annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                           hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                           gaps_col, gaps_row, labels_row, labels_col, ...) 
{
  lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                     ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                     treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                     legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                     annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                     annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                     main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                     fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                     gaps_col = gaps_col, ...)
  res = lo$gt
  mindim = lo$mindim
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = T)
    }
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, 
                                          ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                             res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                          units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                         ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                   ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                    res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
    f(filename, height = height, width = width)
    gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                       border_color = border_color, tree_col = tree_col, 
                       tree_row = tree_row, treeheight_col = treeheight_col, 
                       treeheight_row = treeheight_row, breaks = breaks, 
                       color = color, legend = legend, annotation_col = annotation_col, 
                       annotation_row = annotation_row, annotation_colors = annotation_colors, 
                       annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                       annotation_names_col = annotation_names_col, filename = NA, 
                       main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, hjust_col = hjust_col, 
                       vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                       fontsize_number = fontsize_number, number_color = number_color, 
                       labels_row = labels_row, labels_col = labels_col, 
                       gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    return(gt)
  }
  if (mindim < 3) 
    border_color = NA
  if (!is.na(main)) {
    elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
                          clip = "off")
  }
  if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
    elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
    res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
  }
  if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
    elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }
  elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                fmat, fontsize_number, number_color)
  res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                        name = "matrix")
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                ...)
    elem = do.call(pheatmap:::draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
                          name = "col_names")
  }
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                ...)
    elem = do.call(pheatmap:::draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                          name = "row_names")
  }
  if (!pheatmap:::is.na2(annotation_col)) {
    converted_annotation = convert_annotations(annotation_col, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_col, fontsize, horizontal = T)
    res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
                          name = "col_annotation")
    if (annotation_names_col) {
      elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                              horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                            name = "col_annotation_names")
    }
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    converted_annotation = convert_annotations(annotation_row, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_row, fontsize, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                          name = "row_annotation")
    if (annotation_names_row) {
      elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                              horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                              angle_col = angle_col)
      res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                            name = "row_annotation_names")
    }
  }
  annotation = c(annotation_col[length(annotation_col):1], 
                 annotation_row[length(annotation_row):1])
  annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
  if (length(annotation) > 0 & annotation_legend) {
    elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                             border_color, fontsize = fontsize, ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                          clip = "off", name = "annotation_legend")
  }
  if (!pheatmap:::is.na2(legend)) {
    elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                  ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
                          clip = "off", name = "legend")
  }
  return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                annotation_colors, annotation_legend, annotation_names_row, 
                annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                angle_col, gaps_row, gaps_col, ...) 
{
  if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
    if (!is.null(coln[1])) {
      t = coln
    }
    else {
      t = ""
    }
    tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
    if (annotation_names_row) {
      t = c(t, colnames(annotation_row))
      tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
    }
    longest_coln = which.max(tw)
    gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                fontsize_col, fontsize), ...)
    coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                 rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                  "bigpts")
  }
  else {
    coln_height = unit(5, "bigpts")
  }
  if (!is.null(rown[1])) {
    t = rown
    tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
    if (annotation_names_col) {
      t = c(t, colnames(annotation_col))
      tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
    }
    longest_rown = which.max(tw)
    gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                fontsize_row, fontsize), ...)
    rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                               rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else {
    rown_width = unit(5, "bigpts")
  }
  gp = list(fontsize = fontsize, ...)
  if (!pheatmap:::is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", 
                         textGrob(as.character(names(legend))[longest_break], 
                                  gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                   gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
  }
  else {
    legend_width = unit(0, "bigpts")
  }
  if (is.na(main)) {
    main_height = unit(0, "npc")
  }
  else {
    main_height = unit(1.5, "grobheight", textGrob(main, 
                                                   gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  if (!pheatmap:::is.na2(annotation_col)) {
    annot_col_height = ncol(annotation_col) * (textheight + 
                                                 unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_col_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_col_height = unit(0, "bigpts")
    annot_col_legend_width = unit(0, "bigpts")
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    annot_row_width = ncol(annotation_row) * (textheight + 
                                                unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_row_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_row_width = unit(0, "bigpts")
    annot_row_legend_width = unit(0, "bigpts")
  }
  annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                         "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                         "bigpts")
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - rown_width - legend_width - 
      treeheight_row - annot_row_width - annot_legend_width
  }
  else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
      unit(4, "bigpts")
  }
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - coln_height - 
      treeheight_col - annot_col_height
  }
  else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
      unit(4, "bigpts")
  }
  gt = gtable(widths = unit.c(treeheight_row, rown_width,  
                              mat_width, treeheight_row, legend_width, annot_legend_width), 
              heights = unit.c(main_height, treeheight_col, annot_col_height, 
                               mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
                                                                                    gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                         "bigpts")), "bigpts", valueOnly = T)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/nrow
  mindim = min(cw, ch)
  res = list(gt = gt, mindim = mindim)
  return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
  coord = pheatmap:::find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                 hjust = 1, gp = gpar(...))
  return(res)
}

assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap")









# stacked_barplot <- function(contr.scores, n, data.names, plot.variant = T,
#                             main = NULL, xlab = NULL, ylab = NULL, centering.label = -0.3) {
#   
#   #stopifnot(all(abs(colSums(contr.scores)-1) < 1e-10))
#   
#   k = ncol(contr.scores)
#   
#   data.matrix = matrix(0, nrow = n+1, ncol = k)
#   legend.labels = list()
#   
#   for (ii in 1:k) {
#     
#     contr.lat = top_n_values_indices(contr.scores[,ii], n=n)
#     data.matrix[,ii] = c(contr.lat$values, 1 - sum(contr.lat$values))
#     
#     if (plot.variant) {
#       
#       legend.labels[[ii]] = c(query_locus_names(data.names[contr.lat$indices])$gene, "Rest")
#     }
#     
#     else {
#       legend.labels[[ii]] = c(data.names[contr.lat$indices], "Rest")
#     }
#   }
#   
#   colors = c(colorRampPalette(c("darkblue", "blue", "cyan"))(n), "grey")
#   
#   default.par = par(no.readonly = TRUE)  
#   par(mar = c(1, 4, 2, 2))
#   
#   bar_positions = barplot(data.matrix, beside = FALSE, col = colors, ylim = c(0, 1 + 0.08*n), ylab = "", yaxt = "n", yaxs = "i")
#   
#   for (ii in 1:k) {
#     
#     legend(x = bar_positions[ii] + centering.label, y = 1 + 0.07*n, legend = legend.labels[[ii]], fill = colors, title = paste0("Latent ", ii, " (", round((1-data.matrix[n+1,ii])*100, 1), "%)"))
#   }
#   
#   axis(2, at = seq(0, 1, by = 0.2))
#   mtext(ylab, side = 2, line = 2.5, adj = 0.4 - 0.02*n, cex = 1.4)
#   title(main = main, cex.main = 2.5)
#   
#   par(default.par)
# }






