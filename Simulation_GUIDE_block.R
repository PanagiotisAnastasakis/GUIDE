
source("/Users/panos/Desktop/MSc Thesis/Code/helper.R")
library(ggplot2)
library(pheatmap)
library(grid)


### testing

# pheatmap(W.xl,
#          color = colorRampPalette(c("navy", "gainsboro", "firebrick3"))(100), 
#          #breaks = seq(-max(abs(W.xl.guide)), max(abs(W.xl.guide)), length.out = 101),
#          cluster_rows = F,
#          cluster_cols = F,
#          main = paste0("GUIDE variants-to-latents weights, Large T2D Dataset, L=", k),
#          fontsize = 14)






### Simulation 1: Simple example



set.seed(28)


n = 3000
p = 200
k = 20

W.xl = generate_latent_effect_matrix(n = n, n.latents = k, lower=-2, upper=2, noise.sd = 0.3, dependence.prob = 0.1)
W.lt = t(generate_latent_effect_matrix(n = p, n.latents = k, lower=-2, upper=2, noise.sd = 0.3, dependence.prob = 0.1))


W = W.xl %*% W.lt + matrix(rnorm(n*p, 0, 0.5), nrow = n)


guide.list = get_guide(W, K = k, verbose=F)
W.xl.guide = guide.list$W.xl
W.lt.guide = guide.list$W.lt
weights.guide = rbind(W.xl.guide, W.lt.guide)

weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                            weights.estimated = weights.guide,
                                            n.variants = n)

W.xl.guide.matched = weights.guide.matched$weights.matched.xl
W.lt.guide.matched = weights.guide.matched$weights.matched.lt


degas.list = get_tsvd(scale(W), K = k)
W.xl.degas = degas.list$U
W.lt.degas = degas.list$V
weights.degas = rbind(W.xl.degas, W.lt.degas)

weights.degas.matched = get_matched_weights(weights.true = weights.true,
                                            weights.estimated = weights.degas,
                                            n.variants = n)

W.xl.degas.matched = weights.degas.matched$weights.matched.xl
W.lt.degas.matched = weights.degas.matched$weights.matched.lt





## save at 11.7x9.65


pheatmap(W.xl,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         cluster_rows = F,
         cluster_cols = F,
         main = "",
         fontsize = 15,
         cellheight = 0.21)

grid.text("Simulated Variant Weights", 
          x = 0.5, y = 0.98, gp = gpar(fontsize = 35, fontface = "bold"))

pheatmap(W.xl.guide,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         cluster_rows = F,
         cluster_cols = F,
         main = "",
         fontsize = 15,
         cellheight = 0.21)

grid.text("GUIDE Variant Weights", 
          x = 0.5, y = 0.98, gp = gpar(fontsize = 35, fontface = "bold"))

pheatmap(W.xl.degas,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         cluster_rows = F,
         cluster_cols = F,
         main = "",
         fontsize = 15,
         cellheight = 0.21)

grid.text("DeGAs Variant Weights", 
          x = 0.5, y = 0.98, gp = gpar(fontsize = 35, fontface = "bold"))





## GUIDE

## save at 12.5, 9
plot_simulation(W.xl, W.xl.guide.matched, xlab = "GUIDE Variant Weights", ylab = "Simulated Variant Weights", main = "GUIDE & Simulated Variant Weights")
plot_simulation(t(W.lt), W.lt.guide.matched, xlab = "GUIDE Trait Weights", ylab = "Simulated Trait Weights", main = "GUIDE & Simulated Trait Weights")


#### TSVD

## save at 12.5, 9
plot_simulation(W.xl, W.xl.degas.matched, xlab = "DeGAs Variant Weights", ylab = "Simulated Variant Weights", main = "DeGAs & Simulated Variant Weights")
plot_simulation(t(W.lt), W.lt.degas.matched, xlab = "DeGAs Trait Weights", ylab = "Simulated Trait Weights", main = "DeGAs & Simulated Trait Weights")







#### Simulation 2: arguing for not choosing the GUIDE component selection algorithm

## save them all at 12x9.75

set.seed(22)


### Plot 1



n = 1500
p = 150
n.blocks = 15
K = 50


W.xl = generate_latent_effect_matrix(n = n, n.latents = n.blocks, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1)
W.lt = t(generate_latent_effect_matrix(n = p, n.latents = n.blocks, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1))

W = W.xl %*% W.lt + matrix(rnorm(n*p, 0, 0.5), nrow = n) 


output = get.nlatents(W, starting.K = K, validation.reps = 30, return.estimates = T, alg.typ = "deflation")


ggplot(data.frame(output), aes(x = output)) +
  geom_histogram(binwidth = 0.5, fill = "blue") +
  labs(title = "Histogram of GUIDE Unmixing Matrices Comparison",
       x = "Estimated Components",
       y = "Count" ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min(output), max(output), by = 1)) +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))






### Plots 2,3


n.1 = 200
p.1 = 30

n.2 = 1000
p.2 = 100

n.blocks.1 = 5
n.blocks.2 = 10

K.1 = 10
K.2 = 30

nlat.output.1 = c()
nlat.output.2 = c()

n.reps = 100

for (ii in 1:n.reps) {
  
  
  W.xl = generate_latent_effect_matrix(n = n.1, n.latents = n.blocks.1, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1)
  W.lt = t(generate_latent_effect_matrix(n = p.1, n.latents = n.blocks.1, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1))
  
  W = W.xl %*% W.lt + matrix(rnorm(n.1*p.1, 0, 0.5), nrow = n.1)
  
  K.est.1 = get.nlatents(W, starting.K = K.1, validation.reps = 30, alg.typ = "deflation", verbose = F)$K
  
  nlat.output.1 = c(nlat.output.1, K.est.1)
  
  print(paste("Simulation", ii, "completed."))
}

ggplot(data.frame(nlat.output.1), aes(x = nlat.output.1)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue") +
  labs(title = "GUIDE Component Estimation, Structure 1",
       x = "Estimated Components",
       y = "Count" ) + 
  theme_minimal() +
  scale_x_continuous(breaks = 5:10) + 
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18),
        plot.margin = margin(20, 10, 20, 15))


for (ii in 1:n.reps) {
  
  
  W.xl = generate_latent_effect_matrix(n = n.2, n.latents = n.blocks.2, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1)
  W.lt = t(generate_latent_effect_matrix(n = p.2, n.latents = n.blocks.2, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1))
  
  W = W.xl %*% W.lt + matrix(rnorm(n.2*p.2, 0, 0.5), nrow = n.2) 
  
  K.est.2 = get.nlatents(W, starting.K = K.2, validation.reps = 30, alg.typ = "deflation", verbose = F)$K
  
  nlat.output.2 = c(nlat.output.2, K.est.2)
  
  print(paste("Simulation", ii, "completed."))
}


ggplot(data.frame(nlat.output.2), aes(x = nlat.output.2)) +
  geom_histogram(binwidth = 0.5, fill = "lightblue") +
  labs(title = "GUIDE Component Estimation, Structure 2",
       x = "Estimated Components",
       y = "Count" ) + 
  theme_minimal() +
  scale_x_continuous(breaks = 10:11) +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 18),  
        axis.text.y = element_text(size = 18),
        plot.margin = margin(20, 10, 20, 15))



### Plot 4


W.xl = generate_latent_effect_matrix(n = n.2, n.latents = n.blocks.2, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1)
W.lt = t(generate_latent_effect_matrix(n = p.2, n.latents = n.blocks.2, lower = -2, upper = 2, dependence.prob = 0, noise.sd = 0.1))

W = W.xl %*% W.lt + matrix(rnorm(n.2*p.2, 0, 0.5), nrow = n.2) 


d = get_tsvd(scale(W), K = p.2)$d

ggplot(data.frame(index = 1:length(d), d = d), aes(x = index, y = d)) +
  geom_point(color = "black", size = 2) +       
  labs(title = "Scree plot for SVD on Structure 2",
       x = "Component Number", 
       y = "Singular Value") +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))








## plots 5,6,7,8




n = 1000
p = 100
n.blocks = 10

K = 30

n.reps = 100

dependence_probs = c(0.05, 0.1, 0.15, 0.2)


nlat.outputs = matrix(0, nrow = length(dependence_probs), ncol = n.reps) ## each row is one dependence probability

for (jj in 1:n.reps) {
  
  for (ii in 1:length(dependence_probs)){
    
    W.xl = generate_latent_effect_matrix(n = n, n.latents = n.blocks, lower = -2, upper = 2, dependence.prob = dependence_probs[ii], noise.sd = 0.1)
    W.lt = t(generate_latent_effect_matrix(n = p, n.latents = n.blocks, lower = -2, upper = 2, dependence.prob = dependence_probs[ii], noise.sd = 0.1))
    
    W = W.xl %*% W.lt + matrix(rnorm(n*p, 0, 0.5), nrow = n)
    
    nlat.output = get.nlatents(W, starting.K = K, validation.reps = 30, alg.typ = "deflation", verbose = F)$K
    
    nlat.outputs[ii, jj] = nlat.output
  }
  
  print(paste("Simulation", jj, "completed."))
}



ggplot(data.frame(nlat.outputs[1,]), aes(x = nlat.outputs[1,])) +
  geom_histogram(binwidth = 0.5, fill = "thistle") +
  labs(title = expression(paste("GUIDE Component Estimation, ", italic(P)[dep], " = 0.05")),
       x = "Estimated Components",
       y = "Count" ) +
  scale_x_continuous(breaks = 9:11) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))


ggplot(data.frame(nlat.outputs[2,]), aes(x = nlat.outputs[2,])) +
  geom_histogram(binwidth = 0.5, fill = "thistle") +
  labs(title = expression(paste("GUIDE Component Estimation, ", italic(P)[dep], " = 0.1")),
       x = "Estimated Components",
       y = "Count" ) +
  theme_minimal() +
  scale_x_continuous(breaks = 8:10) +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))

ggplot(data.frame(nlat.outputs[3,]), aes(x = nlat.outputs[3,])) +
  geom_histogram(binwidth = 0.5, fill = "thistle") +
  labs(title = expression(paste("GUIDE Component Estimation, ", italic(P)[dep], " = 0.15")),
       x = "Estimated Components",
       y = "Count" ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))


ggplot(data.frame(nlat.outputs[4,]), aes(x = nlat.outputs[4,])) +
  geom_histogram(binwidth = 0.5, fill = "thistle") +
  labs(title = expression(paste("GUIDE Component Estimation, ", italic(P)[dep], " = 0.2")),
       x = "Estimated Components",
       y = "Count" ) +
  theme_minimal() +
  scale_x_continuous(breaks = 2:7) +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))






