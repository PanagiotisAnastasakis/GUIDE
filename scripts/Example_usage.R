

## To demonstrate how GUIDE and DeGAs works, we generate data from a polygenic 
## additive model for 1000 variants, 100 traits and 15 latent components, where
## the weights have low sparsity and are sampled from the Laplace distribution.

n = 1000
p = 100
K = 15

sim_polygenic = polygenic_sim(n = n, p = p, k = K,
                              sampler = "laplace",
                              sparsity = "low")

W = sim_polygenic$W
W.xl = sim_polygenic$W.xl
W.lt = sim_polygenic$W.lt

weights.true = rbind(W.xl, W.lt)

## To obtain the GUIDE and DeGAs decompositions, we use the functions 'get_guide'
## and 'get_degas' respectively.


## For GUIDE, the function 'get_guide' incorporates the 
## proposed ICASSO extension, where multiple ICA runs are considered 
## to get the final estimate. To use ICASSO with a specified number of runs
## (e.g. 50 below), one must set the 'ica_runs' argument accordingly. If it is set to 1,
## then the standard GUIDE method is used, where ICA is run once.
## Note that all data preprocessing for GUIDE (mean-centering and whitening) is performed automatically.

guide.dec = get_guide(B = W, K = K, ica_runs = 50) ## 'B' is the summary statistics matrix and 'K' the number of components

## Here, 'guide.dec' is a list of outputs related to the implementation of GUIDE,
## described in the file GUIDE.R.


W.xl.guide = guide.dec$W.xl ## the variant-to-latent weights (components are in the columns)
W.lt.guide = guide.dec$W.lt ## the latent-to-trait weights (components are in the columns)

weights.guide = rbind(W.xl.guide, W.lt.guide)


## Using the function 'get_matched_weights' (assuming that the true weights exist and are known),
## one can change the ordering of the simulated weights to match the true ones, also
## correcting for possible sign changes.

weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                            weights.estimated = weights.guide,
                                            n.variants = n)

W.xl.guide.matched = weights.guide.matched$weights.matched.xl
W.lt.guide.matched = weights.guide.matched$weights.matched.lt


## Using the function 'plot_simulation', we can plot the estimated weights against
## the true ones for each component:

plot_simulation(W.xl, W.xl.guide.matched, xlab = "GUIDE Variant Weights", ylab = "Simulated Variant Weights", main = "GUIDE & Simulated Variant Weights")
plot_simulation(W.lt, W.lt.guide.matched, xlab = "GUIDE Trait Weights", ylab = "Simulated Trait Weights", main = "GUIDE & Simulated Trait Weights")


## Using the output of 'get_guide', we can also compute the Cluster Quality Index (CQI)
##scores of each component note that this is only available with the ICASSO extension,
## meaning when 'ica_runs' is greater than 1. These scores could then be used to sort
## the components according to their significance.

get_cqi_values(guide.dec$cors.unmix, guide.dec$clusters)


## we can further explore the clustering of ICASSO within GUIDE by plotting the 
## dendrogram of the clustering process


plot(guide.dec$hc,
     hang = 0.01,    
     labels = F,
     main = "Dendrogram of ICASSO Clustering",
     ylab = "Distance",
     xlab = "", sub = "",
     las=1, cex.axis = 1.7, cex.lab=3, cex.main = 3)


## Additional analyses can be performed with the resulting GUIDE weights W.xl.guide,
## W.lt.guide, for instance computing the kurtosis of each component, visualizing the 
## most significant elements of a component (variants or traits), as determined by
## the magnitude of their weights, or comparing the overall quality of the decomposition
## with that of other methods, such as DeGAs.

## To run DeGAs, which applies SVD on a given summary statistics matrix, we can 
## use the function 'get_tsvd'.


degas.list = get_tsvd(B = scale(W), K = K) ## the arguments are as above

W.xl.degas = degas.list$U ## the variant-to-latent weights
W.lt.degas = degas.list$V ## the latent-to-trait weights

## with DeGAs, we can use the singular values to assess component significance

W.lt.degas = degas.list$d


## Unlike GUIDE, DeGAs supports methods for choosing a number of components using
## popular statistical approaches related to SVD. In this implementation, two such
## approaches are supported: the above average eigenvalue criterion and a specific percentage
## of variance to be explained by the components.

## To use the above average eigenvalue criterion, we must have K=0 and set the argument
## 'method' to 'average_eig':

degas.list.average_eig = get_tsvd(B = scale(W), K = 0, method = "average_eig")

## To choose components according to a specific percentage of variance explained,
## again K must be set to 0 and now the argument 'method' is set to 'var_explained',
## while we must also specify the percentage 'p' using a value from 0 to 1.

degas.list.var_explained = get_tsvd(B = scale(W), K = 0, method = "var_explained", p = 0.8)


## The output is the same as before, the only difference being the statistical estimation of 
## the number of components, instead of the use of a fixed value.

