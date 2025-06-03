
library(ExtDist)
library(ggplot2)


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





### Function to generate block latent components with a given level of dependency

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
      remaining.effs = replicate(n - n/n.latents, ifelse(runif(1) < dependence.prob, runif(1, min = lower, max = upper), 0))
                                 
      #                            {
      #   
      #   if (runif(1) < dependence.prob) {
      #     runif(1, min = lower, max = upper)
      #   } else {
      #     0
      #   }
      # })
      
      W[-(row_start:row_end), i] = remaining.effs
    }
  }
  
  W = W + matrix(rnorm(n*n.latents, 0, noise.sd), nrow = n, ncol = n.latents) ## add noise
  
  return(W)
}



## Function for simulations under a polygenic additive model


polygenic_sim <- function(n, p, k, sampler = "gaussian", sparsity = "none", h.sq = 0.1, Nsamples = 1e5) {#, return.afs = F) {
  
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
  
  #if (return.afs) return(list(W = W, W.xl = W.xl, W.lt = W.lt, afs = afs))
  
  return(list(W = W, W.xl = W.xl, W.lt = W.lt, afs = afs))
}

