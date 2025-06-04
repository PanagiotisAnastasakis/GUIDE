
## Function to get a GWAS summary statistics matrix decomposition based on the method
## Decomposition of Genetic Associations (DeGAs), which employs truncated SVD.
## Two approaches for choosing the number of components are considered:

## 1) Proportion of variance explained: we require that the first K components explain p*100% of the variance
## 2) Above average eigenvalues. The "average" is computed as the sum of squared singular values, divided by the rank of the matrix


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

