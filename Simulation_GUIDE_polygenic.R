
library(ggplot2)
library(pheatmap)
library(vegan)
library(matrixcalc)
library(moments)
library(ExtDist)
library(patchwork)
library(latex2exp)
library(cowplot)




### Simulation 1

### We simulate z-score matrices consisting of different number of variants and traits,
### using the Gaussian and Laplace Sampler and we consider 3 different levels of
### sparsity; None, Low, High.



set.seed(4224)


## (SNPs, Traits, Latents)

#conf.1 = c(100, 20, 4)
conf.1 = c(500, 50, 5)
conf.2 = c(1000, 100, 10)
#conf.4 = c(3000, 150, 20)


n.reps = 300

confs = list(conf.1, conf.2)


results = list("C1" = list(), "C2" = list())

for (i in 1:length(results)) {
  
  configuration = names(results)[i]
  results[[configuration]] = list("gaussian" = list(), "laplace" = list())

  for (sampler in names(results[[configuration]])) {
    results[[configuration]][[sampler]] = list("none" = list(), "low" = list(), "high" = list())
    
    for (sparsity in names(results[[configuration]][[sampler]])) {
      
      avg.cor.guide.xl = c()
      avg.cor.guide.lt = c()
      
      avg.cor.degas.xl = c()
      avg.cor.degas.lt = c()
      
      avg.kurtosis.guide = c()
      avg.kurtosis.degas = c()
      
      n = confs[[i]][1]
      p = confs[[i]][2]
      k = confs[[i]][3]
      
      for (jj in 1:n.reps) {
        
        sim_polygenic = polygenic_sim(n = n, p = p, k = k, sampler = sampler, sparsity = sparsity)
        
        W = sim_polygenic$W
        W.xl = sim_polygenic$W.xl
        W.lt = sim_polygenic$W.lt
        
        weights.true = rbind(W.xl, W.lt)
        
        ## GUIDE
        
        #ica_clustering = guide_icasso(W = W, K = k, reps = 20)
        #a = ica_clustering$optimal.unmixing.matrix
        guide.list = get_guide(W, K = k, ica_runs = 20, verbose = F)
        
        W.xl.guide = guide.list$W.xl
        W.lt.guide = guide.list$W.lt
        weights.guide = rbind(W.xl.guide, W.lt.guide)
        
        weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                                    weights.estimated = weights.guide,
                                                    n.variants = n)
        
        W.xl.guide.matched = weights.guide.matched$weights.matched.xl
        W.lt.guide.matched = weights.guide.matched$weights.matched.lt
        
        
        ## DeGAs
        
        degas.list = get_tsvd(scale(W), K = k)
        W.xl.degas = degas.list$U
        W.lt.degas = degas.list$V
        weights.degas = rbind(W.xl.degas, W.lt.degas)
        
        weights.degas.matched = get_matched_weights(weights.true = weights.true,
                                                    weights.estimated = weights.degas,
                                                    n.variants = n)
        
        W.xl.degas.matched = weights.degas.matched$weights.matched.xl
        W.lt.degas.matched = weights.degas.matched$weights.matched.lt
        
        
        ## Statistics computed for both
        
        avg.cor.guide.xl = c(avg.cor.guide.xl, mean(sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.guide.matched[,ii]))))
        avg.cor.guide.lt = c(avg.cor.guide.lt, mean(sapply(1:k, function(ii) cor(W.lt[,ii], W.lt.guide.matched[,ii]))))
        
        avg.cor.degas.xl = c(avg.cor.degas.xl, mean(sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.degas.matched[,ii]))))
        avg.cor.degas.lt = c(avg.cor.degas.lt, mean(sapply(1:k, function(ii) cor(W.lt[,ii], W.lt.degas.matched[,ii]))))
        
        avg.kurtosis.guide = c(avg.kurtosis.guide, mean(apply(weights.guide, 2, kurtosis)) - 3)
        avg.kurtosis.degas = c(avg.kurtosis.degas, mean(apply(weights.degas, 2, kurtosis)) - 3)
        
      }
      
      resulting.stats.guide = data.frame(avg.cor.xl = avg.cor.guide.xl,
                                         avg.cor.lt = avg.cor.guide.lt,
                                         avg.kurtosis = avg.kurtosis.guide)
      
      resulting.stats.degas = data.frame(avg.cor.xl = avg.cor.degas.xl,
                                         avg.cor.lt = avg.cor.degas.lt,
                                         avg.kurtosis = avg.kurtosis.degas)
      
      
      results[[configuration]][[sampler]][[sparsity]][["GUIDE"]] = resulting.stats.guide
      results[[configuration]][[sampler]][[sparsity]][["DeGAs"]] = resulting.stats.degas
      
      
      print(paste("Simulation for configuration", configuration, "with", sampler, "sampler and", sparsity, "sparsity completed."))
    }
  }
}



saveRDS(results, file = "polygenic_results.rds")




### Things to write about

### How can TSVD be viewed as a preprocessing step for GUIDE (whitening). Show it mathematically. 
### Discuss the results in (some) detail





### Simulation 2

### We evaluate the performance of GUIDE and DeGAs under Model Misspecification

## We compute the average correlations by keeping only the min(k_true, k_fitted) components and averaging, same as before.

## if GUIDE is still able to find the true components, we can argue that in a practical setting,
## in absence of a component selection strategy, it would be preferable to select a higher number
## as then the meaningful components can still be found.

### We consider one configuration consisting of 1000 variants, 100 traits and 10 components.
### We then compute GUIDE and DeGAs for k=2,3,...,30, thus evaluating both over-specification
### and under-specification

### We test both Gaussian & Laplace Sampler, across low and high sparsity.

set.seed(2)

n.reps = 200

n = 1000
p = 100
k = 10

k.est.vals = 2:30


results_misspec = list("gaussian" = list(), "laplace" = list())
sim_list = list("gaussian" = list(), "laplace" = list())
  
for (sampler in names(results_misspec)) {

  results_misspec[[sampler]] = list("low" = list(), "high" = list())
  sim_list[[sampler]] = list("low" = list(), "high" = list())
  
  for (sparsity in names(results_misspec[[sampler]])) {
    
    results_misspec[[sampler]][[sparsity]] = vector("list", length(k.est.vals)+1) ## the first list corresponds to k=1 and is empty
    
    sim_list[[sampler]][[sparsity]] = lapply(1:n.reps, function(...) polygenic_sim(n = n, p = p, k = k, sampler = sampler, sparsity = sparsity)) ## using the same simulated matrices
    
    for (k.est in k.est.vals) {
      
      avg.cor.guide.xl = c()
      avg.cor.guide.lt = c()
      avg.cor.guide.total = c()
      cor.min.guide = c()

      avg.cor.degas.xl = c()
      avg.cor.degas.lt = c()
      avg.cor.degas.total = c()
      
      for (jj in 1:n.reps) {
          
        sim_polygenic = sim_list[[sampler]][[sparsity]][[jj]]
        
        W = sim_polygenic$W
        W.xl = sim_polygenic$W.xl
        W.lt = sim_polygenic$W.lt
        
        weights.true = rbind(W.xl, W.lt)
        
        ## GUIDE
        #ica_clustering = guide_icasso(W = W, K = k.est, reps = 20)
        #a =ica_clustering$optimal.unmixing.matrix
        guide.list = get_guide(W, K = k.est, ica_runs = 20, verbose = F)
        
        W.xl.guide = guide.list$W.xl
        W.lt.guide = guide.list$W.lt
        weights.guide = rbind(W.xl.guide, W.lt.guide)
        
        ## DeGAs
        
        degas.list = get_tsvd(scale(W), K = k.est)
        W.xl.degas = degas.list$U
        W.lt.degas = degas.list$V
        weights.degas = rbind(W.xl.degas, W.lt.degas)
        
        if (k.est < k) {
          
          weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                                      weights.estimated = weights.guide,
                                                      n.variants = n)
          
          W.xl.matched = weights.guide.matched$weights.matched.xl
          W.lt.matched = weights.guide.matched$weights.matched.lt
          
          weights.true.matched = rbind(W.xl.matched, W.lt.matched)
          
          avg.cor.guide.xl = c(avg.cor.guide.xl, mean(sapply(1:k.est, function(ii) cor(W.xl.matched[,ii], W.xl.guide[,ii]))))
          avg.cor.guide.lt = c(avg.cor.guide.lt, mean(sapply(1:k.est, function(ii) cor(W.lt.matched[,ii], W.lt.guide[,ii]))))
          
          matched.cors = sapply(1:k.est, function(ii) cor(weights.true.matched[,ii], weights.guide[,ii]))
          cor.min.guide = c(cor.min.guide, min(matched.cors))
          avg.cor.guide.total = c(avg.cor.guide.total, mean(matched.cors))
          
          
          weights.degas.matched = get_matched_weights(weights.true = weights.true,
                                                      weights.estimated = weights.degas,
                                                      n.variants = n)
          
          W.xl.matched = weights.degas.matched$weights.matched.xl
          W.lt.matched = weights.degas.matched$weights.matched.lt
          
          weights.true.matched = rbind(W.xl.matched, W.lt.matched)
          
          matched.cors = sapply(1:k.est, function(ii) cor(weights.true.matched[,ii], weights.degas[,ii]))
          avg.cor.degas.total = c(avg.cor.degas.total, mean(matched.cors))
          
          avg.cor.degas.xl = c(avg.cor.degas.xl, mean(sapply(1:k.est, function(ii) cor(W.xl.matched[,ii], W.xl.degas[,ii]))))
          avg.cor.degas.lt = c(avg.cor.degas.lt, mean(sapply(1:k.est, function(ii) cor(W.lt.matched[,ii], W.lt.degas[,ii]))))
          
        }
          
        else {
            
          weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                                      weights.estimated = weights.guide,
                                                      n.variants = n)
          
          W.xl.guide.matched = weights.guide.matched$weights.matched.xl
          W.lt.guide.matched = weights.guide.matched$weights.matched.lt
          
          weights.guide.matched.total = rbind(W.xl.guide.matched, W.lt.guide.matched)
          
          avg.cor.guide.xl = c(avg.cor.guide.xl, mean(sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.guide.matched[,ii]))))
          avg.cor.guide.lt = c(avg.cor.guide.lt, mean(sapply(1:k, function(ii) cor(W.lt[,ii], W.lt.guide.matched[,ii]))))
          
          matched.cors = sapply(1:k, function(ii) cor(weights.true[,ii], weights.guide.matched.total[,ii]))
          cor.min.guide = c(cor.min.guide, min(matched.cors))
          avg.cor.guide.total = c(avg.cor.guide.total, mean(matched.cors))
          
          
          weights.degas.matched = get_matched_weights(weights.true = weights.true,
                                                      weights.estimated = weights.degas,
                                                      n.variants = n)
          
          W.xl.degas.matched = weights.degas.matched$weights.matched.xl
          W.lt.degas.matched = weights.degas.matched$weights.matched.lt
          
          weights.degas.matched.total = rbind(W.xl.degas.matched, W.lt.degas.matched)
          matched.cors = sapply(1:k, function(ii) cor(weights.true[,ii], weights.degas.matched.total[,ii]))
          
          avg.cor.degas.total = c(avg.cor.degas.total, mean(matched.cors))
          
          avg.cor.degas.xl = c(avg.cor.degas.xl, mean(sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.degas.matched[,ii]))))
          avg.cor.degas.lt = c(avg.cor.degas.lt, mean(sapply(1:k, function(ii) cor(W.lt[,ii], W.lt.degas.matched[,ii]))))
          
        }
      }
      
    results_misspec[[sampler]][[sparsity]][[k.est]][["GUIDE"]] = data.frame(avg.cor.xl = avg.cor.guide.xl,
                                                                            avg.cor.lt = avg.cor.guide.lt,
                                                                            avg.cor.total = avg.cor.guide.total,
                                                                            cor.min = cor.min.guide)
    
    results_misspec[[sampler]][[sparsity]][[k.est]][["DeGAs"]] = data.frame(avg.cor.xl = avg.cor.degas.xl,
                                                                            avg.cor.lt = avg.cor.degas.lt,
                                                                            avg.cor.total = avg.cor.degas.total)

    print(paste("Simulation for k =", k.est, "with", sampler, "sampler and", sparsity, "sparsity completed."))
    }
  }
}


saveRDS(results_misspec, file = "polygenic_results_misspec.rds")








### Simulation 3

## We can consider the same number of variants, traits anc components as before, and
## we overestimate the true number of components by 10.


set.seed(20)

n = 1000
p = 100
k = 10
k.est = 20



n.reps = 200


results_misspec_recovery_guide = list("gaussian" = list(), "laplace" = list())

for (sampler in names(results_misspec_recovery_guide)) {
  
  results_misspec_recovery_guide[[sampler]] = list("low" = list(), "high" = list())
  
  for (sparsity in names(results_misspec_recovery_guide[[sampler]])) {
    
    n.captured.components = matrix(NA, nrow = n.reps, ncol = k.est)
    n.captured.components.top5 = matrix(NA, nrow = n.reps, ncol = k.est)
    
    for (jj in 1:n.reps) {
      
      sim_polygenic = polygenic_sim(n = n, p = p, k = k, sampler = sampler, sparsity = sparsity)
      
      #sim_polygenic = sim_list[[sampler]][[sparsity]][[jj]] ## using the generated data from simulation 2
      
      W = sim_polygenic$W
      W.xl = sim_polygenic$W.xl
      W.lt = sim_polygenic$W.lt
      
      weights.true = rbind(W.xl, W.lt)
      
      #ica_clustering = guide_icasso(W = W, K = k.est, reps = 20)
      #a = ica_clustering$optimal.unmixing.matrix
      guide.list = get_guide(W, K = k.est, ica_runs = 20, verbose = F)
      
      W.xl.guide = guide.list$W.xl
      W.lt.guide = guide.list$W.lt
      weights.guide = rbind(W.xl.guide, W.lt.guide)
      cqi.vals = get_cqi_values(guide.list$cors.unmix, guide.list$clusters)
      
      matches = get_matched_weights(weights.true, weights.guide, n.variants = n,
                                      return.only.indices = TRUE)[,2] ## the matched components
      
      n.captured.components[jj,] = sapply(1:k.est, function(k.est.ord) {
        length(intersect(matches, order(cqi.vals, decreasing = TRUE)[1:k.est.ord]))
      })
      
      n.captured.components.top5[jj,] = sapply(1:k.est, function(k.est.ord) {
        length(intersect(matches[1:5], order(cqi.vals, decreasing = TRUE)[1:k.est.ord]))
      })
      
      
    }
      
    results_misspec_recovery_guide[[sampler]][[sparsity]][["all"]] = n.captured.components
    results_misspec_recovery_guide[[sampler]][[sparsity]][["top5"]] = n.captured.components.top5
    
    
    print(paste("Simulation for", sampler, "sampler and", sparsity, "sparsity completed."))
  
  }
}


saveRDS(results_misspec_recovery_guide, file = "results_misspec_recovery_guide.rds")



