



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