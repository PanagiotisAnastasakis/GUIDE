
library(ExtDist)
library(ggplot2)


## Function to find the optimal permutation (and possible sign change) to match
## a set of true component weights and a set of estimated ones by maximizing their correlations.

get_matched_weights <- function(weights.true, weights.estimated, n.variants,
                                return.only.indices = FALSE) {
  
  ### Inputs
  ## weights.true: the simulated (true) variant and trait weights to the latent components concatenated. It has dimension (n+p)xL,
  ##               where n is the number of variants, p is the number of traits and L is the true number of components.
  ## weights.estimated: the estimated variant and trait weights to the latent factors concatenated, with dimension (n+p)xK,
  ##               where K is the number of components they are estimated with and may differ from L.
  ## n.variants: the number of variants in the genetic architecture (used to separate the variant and trait weights in the end)
  ## return.only.indices: whether to return only the matched indices between the two sets of weights, and not the weights themselves
  
  ### Output:
  ## If return.only.indices == TRUE, then only the indices for the matchings are returned
  ## If L>K, then the first K true weights that are best matched with all the K estimated ones are returned as a list, separately for variants and traits
  ## If L<=K, then all the L true weights are matched with L of the estimated ones, which are then returned as a list, separately for variants and traits
  
  L = ncol(weights.true) 
  K = ncol(weights.estimated) 
  
  ## Computing the correlation between each set of true and estimated component weights:
  
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
    
    weights.matched = weights.true[,sorted.matches[,1]] %*% diag(sorted.matches[,3]) ## we re-order the estimated components and multiply them with either -1 or 1 to correct for sign change
    
    weights.matched.xl = head(weights.matched, n.variants) ## If there are fewer estimated weights than true ones, then we match a subset of the true components with all estimated ones!! 
                                                           ## Also, we keep the first n rows, as these correspond to variants (W.xl.estimated) and the rest (below) are the traits (W.lt.estimated)
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





## Function to plot either the variants-to-latents or the latents-to-traits weights, both assumed to have been matched
## so that columns with the same index correspond to the same component.

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





### Function to generate summary statistics of block genetic architectures with a given level of dependency

block_sim <- function(n, p, k, lower, upper, dependence.prob = 0.1,
                                          noise.sd.1 = 0.02, noise.sd.2 = 0.02) {
  
  ### Inputs
  ## n: the number of vertices (either variants or traits)
  ## n.latents: the number of latent components of the simulated structure
  ## lower, upper: the lower and upper bounds for the Uniform distribution from which the weights are sampled
  ## dependence.prob: value between 0 and 1, corresponding to the degree of dependence between components
  ## noise.sd.1: the standard deviation for noise from the gaussian distribution, added to the variant and trait weights
  ## noise.sd.1: the standard deviation for noise from the gaussian distribution, added to the final block matrix
  
  ### Output: a list consisting of the simulated summary statistics matrix and the variant and trait weights
  
  W.xl = matrix(0, nrow = n, ncol = k)
  W.lt = matrix(0, nrow = p, ncol = k)
  
  for (i in seq_len(k)) {
    
    row_start_xl = (i - 1) * (n / k) + 1
    row_end_xl = i * (n / k)
    block_xl = runif(n / k, min = lower, max = upper)
    W.xl[row_start_xl:row_end_xl, i] = block_xl
    
    row_start_lt = (i - 1) * (p / k) + 1
    row_end_lt = i * (p / k)
    block_lt = runif(p / k, min = lower, max = upper)
    W.lt[row_start_lt:row_end_lt, i] = block_lt
    
    if (dependence.prob > 0) {
      
      remaining.effs.xl = replicate(n - n/k, ifelse(runif(1) < dependence.prob, runif(1, min = lower, max = upper), 0))
      W.xl[-(row_start_xl:row_end_xl), i] = remaining.effs.xl
      
      remaining.effs.lt = replicate(p - p/k, ifelse(runif(1) < dependence.prob, runif(1, min = lower, max = upper), 0))
      W.lt[-(row_start_lt:row_end_lt), i] = remaining.effs.lt
    }
  }
  
  W.xl = W.xl + matrix(rnorm(n*k, 0, noise.sd.1), nrow = n, ncol = k) 
  W.lt = W.lt + matrix(rnorm(p*k, 0, noise.sd.1), nrow = p, ncol = k) 
  
  W = W.xl %*% t(W.lt) + matrix(rnorm(n*p, 0, noise.sd.2), nrow = n)
  
  return(list(W = W, W.xl = W.xl, W.lt = W.lt))
}



## Function for simulating summary statistics under a polygenic additive model


polygenic_sim <- function(n, p, k, sampler = "gaussian", sparsity = "none", h.sq = 0.1, Nsamples = 1e5) {
  
  ### Inputs
  ## n,p,k: the number of variants, traits, components respectively
  ## sampler: the distribution from which the weights are sampled from (either gaussian or laplace)
  ## sparsity: the level of sparsity for the weights, which is none, low, or high
  ## h.sq: the heritability assumed in the genetic architecture
  ## Nsamples: the number of samples assumed for each effect size. Here it controls the power and thus the standard error. Higher sample size means smaller SE
  
  ### Output: a list consisting of the simulated summary statistics matrix, the variant and trait weights, as well as the simulated allele frequencies for each variant
  
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
  
  E = do.call("cbind", lapply(1:p, function(j) rnorm(n, sd = SEs))) ## noise (assuming independent SNPs, hence LD=0)
  
  B = B + E ## adding noise
  
  W = B / replicate(p, SEs) ## dividing by the standard errors to get the z-scores
  

  return(list(W = W, W.xl = W.xl, W.lt = W.lt, afs = afs))
}

