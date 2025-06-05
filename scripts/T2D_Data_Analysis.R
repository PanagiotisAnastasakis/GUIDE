
## This file contains the script for analyzing T2D data from the paper 
## "Multiancestry polygenic mechanisms of type 2 diabetes", by K. Smith et al.


library(ggplot2)
library(cowplot)
library(pheatmap)
library(moments)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggbeeswarm)


set.seed(6941125)


K = 12 ## Number of clusters used in the bNMF analysis of this data. We use the same for both GUIDE and DeGAs

def_par = par(no.readonly = TRUE)


## Running DeGAs

degas.list = get_tsvd(scale(df), K = K)

W.xl.degas = degas.list$U
W.lt.degas = degas.list$V

W.xl.degas_signif = W.xl.degas
W.xl.degas_signif[W.xl.degas_signif^2 < 1/nrow(df)] = 0

W.lt.degas_signif = W.lt.degas
W.lt.degas_signif[W.lt.degas_signif^2 < 1/ncol(df)] = 0


## Running GUIDE with ICASSO

ica_runs = 100

guide.list = get_guide(df, K = K, ica_runs = ica_runs)

cqi_values = get_cqi_values(guide.list$cors.unmix, guide.list$clusters)

cluster_ordering = order(cqi_values, decreasing = TRUE)

W.xl.guide = guide.list$W.xl
W.lt.guide = guide.list$W.lt

W.xl.guide = W.xl.guide[,cluster_ordering]
W.lt.guide = W.lt.guide[,cluster_ordering]

W.xl.guide_signif = W.xl.guide
W.xl.guide_signif[W.xl.guide_signif^2 < 1/nrow(df)] = 0

W.lt.guide_signif = W.lt.guide
W.lt.guide_signif[W.lt.guide_signif^2 < 1/ncol(df)] = 0




## plotting clustering results

## CQI scores


data = data.frame(x = 1:K, score = sort(cqi_values, decreasing = T))

## save at 19.5x11

data %>%
  ggplot(aes(x = x, y = score)) +
  geom_line(color = "grey") +
  geom_point(shape = 21, color = "black", fill = "steelblue", size = 6) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +  # Remove minor x grid lines
  scale_x_continuous(breaks = 1:K) +
  scale_y_continuous(
    breaks = seq(0.5, 1, by = 0.1),
    limits = c(min(cqi_values) - 0.05, NA)
  )+
  labs(title = "CQI Scores Across GUIDE Components",
          x = "Component Index",
          y = "CQI") +
  theme(plot.title = element_text(size = 30, hjust = 0.5, vjust = 3),
      axis.title.x = element_text(size = 24, vjust = -3),
      axis.title.y = element_text(size = 24, vjust = 3),
      axis.text.x = element_text(size = 16),  
      axis.text.y = element_text(size = 16),
      plot.margin = margin(20, 15, 25, 15))


## dendrogram

cut_height = sort(guide.list$hc$height, decreasing = TRUE)[K-1] 

par(mar = c(1, 7, 4, 0), mgp = c(4.5, 1, 0))  

## save at 16x11
plot(guide.list$hc,
     hang = 0.01,    
     labels = F,
     main = "Dendrogram of ICASSO Clustering",
     ylab = "Distance",
     xlab = "", sub = "",
     las=1, cex.axis = 1.7, cex.lab=3, cex.main = 3)

segments(x0 = 0.5, 
         x1 = K*ica_nruns + 10, 
         y0 = cut_height - 0.005, 
         y1 = cut_height - 0.005, 
         col = "red", 
         lwd = 2, 
         lty = 2)

par(def_par)




## showing that one ICA run is less consistent than the proposed ICASSO extension

n_runs = 150

d = get_nlatents(df, starting.K = K, validation.reps = n_runs, return.estimates = T)

sum(d == K)/(n_runs*(n_runs-1)/2)

## save at 13x10

single_ICA_comps = ggplot(data.frame(d), aes(x = d)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 0.5, fill = "lightcoral") +
  labs(title = "Single ICA Runs",
       x = "Matched Components ",
       y = "Frequency" ) +
  theme_minimal() +
  scale_x_continuous(breaks = min(d):max(d)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 17, vjust = -0.5),              
        axis.title.y = element_text(size = 17, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))



### ICASSO

unmixing.matrices = list()

for (ii in 1:n_runs) {

  unmix_matrix_guide = get_guide(df, K = K, ica_runs = 100)$A
  unmixing.matrices[[ii]] = unmix_matrix_guide
  print(ii)
}


combs = combn(1:n_runs, 2)

n.selected.components = c()

for (i in 1:ncol(combs)) {
  
  idx1 = combs[1, i]
  idx2 = combs[2, i]
  
  A.1 = unmixing.matrices[[idx1]]
  A.2 = unmixing.matrices[[idx2]]
  
  components = 0
  
  for (ii in 1:K) {
    for (jj in 1:K) {
      
      if (abs(cor(A.1[,ii], A.2[,jj])) > 0.95) { 
        
        components = components + 1
        break
      }
    }
  }
  n.selected.components = c(n.selected.components, components)
}




icasso_comps = ggplot(data.frame(n.selected.components), aes(x = n.selected.components)) +
  geom_histogram(aes(y = after_stat(count / sum(count))),
                 binwidth = 0.5, fill = "lightcoral") +
  labs(title = "ICASSO",
       x = "Matched Components ",
       y = "Frequency" ) +
  theme_minimal() +
  scale_x_continuous(breaks = 3:12, limits = c(3, 12.25)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 17, vjust = -0.5),              
        axis.title.y = element_text(size = 17, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15))




plots_icasso_ica <- plot_grid(
  single_ICA_comps, icasso_comps,
  ncol = 2, align = "hv"
)

icasso_ica_total <- ggdraw() +
  draw_label("Matched Components per Unmixing Matrix Comparison", fontface = 'plain', size = 28, hjust = 0.5, vjust = -10.6) +
  draw_plot(
    plot_grid(
      plots_icasso_ica + theme(plot.margin = margin(10, -90, 25, 15)),
      ncol = 2, rel_widths = c(1, 0.12)
    ),
    y = -0.08
  )

## save at 16x6.9 inches
icasso_ica_total







#### Cluster plots

## Circle plots 

## Getting variant and gene names


variant_names = rownames(df)
gene_names = query_locus_names(variant_names)$gene

phenotype_names = colnames(df)

phenotype_names = gsub("_female$", "_F", phenotype_names, ignore.case = TRUE)
phenotype_names = gsub("_male$", "_M", phenotype_names, ignore.case = TRUE)
phenotype_names = gsub(".female$", "_F", phenotype_names, ignore.case = TRUE)
phenotype_names = gsub(".male$", "_M", phenotype_names, ignore.case = TRUE)

phenotype_names[phenotype_names == "Haemoglobin.concentration_M"] = "Hgb.concentration_M"
phenotype_names[phenotype_names == "Haemoglobin.concentration_F"] = "Hgb.concentration_F"
phenotype_names[phenotype_names == "High.light.scatter.reticulocyte.count_M"] = "HL_Retics_M"
phenotype_names[phenotype_names == "High.light.scatter.reticulocyte.count_F"] = "HL Retics_F"

total_names = c(gene_names, phenotype_names)

total_weights_guide = rbind(W.xl.guide_signif, W.lt.guide_signif) 


cluster_number = 12

## save at 13x13
cluster_circle_plot(total_weights_guide[,cluster_number],
                    total_names = total_names,
                    title = paste("GUIDE Cluster", cluster_number),
                    rotate = 20)
## 6 rotate 60

## 1 -> cholesterol
## 2 -> Beta Cell 2
## 3 -> Beta Cell 1
## 4 -> Bilirubin
## 5 -> ??
## 6 -> ALP Negative
## 7 -> Liver-Lipid
## 8 -> Lipodystrophy 2
## 9 -> Proinsulin (loosely, ony some traits)
## 10 -> ??
## 11 -> Obesity
## 12 -> Lipodystrophy 1




total_weights_degas = rbind(W.xl.degas_signif, W.lt.degas_signif) 

cluster_number = 12

## save at 13x13
cluster_circle_plot(total_weights_guide[,1],
                    total_names = total_names,
                    title = paste("GUIDE Cluster", 1),
                    rotate = 20)






## Jaccard Indices


#### Jaccard matrix heatmap of GUIDE against DeGAs clusters


jaccard_guide_degas_xl = compute_jaccard_matrix(abs(W.xl.guide_signif), abs(W.xl.degas_signif))
jaccard_guide_degas_lt = compute_jaccard_matrix(abs(W.lt.guide_signif), abs(W.lt.degas_signif))

jaccard_guide_xl = compute_jaccard_matrix(abs(W.xl.guide_signif), abs(W.xl.guide_signif))
jaccard_guide_lt = compute_jaccard_matrix(abs(W.lt.guide_signif), abs(W.lt.guide_signif))

jaccard_degas_xl = compute_jaccard_matrix(abs(W.xl.degas_signif), abs(W.xl.degas_signif))
jaccard_degas_lt = compute_jaccard_matrix(abs(W.lt.degas_signif), abs(W.lt.degas_signif))


rownames(jaccard_guide_degas_xl) = paste0("GUIDE_", 1:K)
colnames(jaccard_guide_degas_xl) = paste0("DeGAs_", 1:K)

rownames(jaccard_guide_degas_lt) = paste0("GUIDE_", 1:K)
colnames(jaccard_guide_degas_lt) = paste0("DeGAs_", 1:K)


rownames(jaccard_guide_xl) = paste0("GUIDE_", 1:K)
colnames(jaccard_guide_xl) = paste0("GUIDE_", 1:K)

rownames(jaccard_guide_lt) = paste0("GUIDE_", 1:K)
colnames(jaccard_guide_lt) = paste0("GUIDE_", 1:K)


rownames(jaccard_degas_xl) = paste0("DeGAs_", 1:K)
colnames(jaccard_degas_xl) = paste0("DeGAs_", 1:K)

rownames(jaccard_degas_lt) = paste0("DeGAs_", 1:K)
colnames(jaccard_degas_lt) = paste0("DeGAs_", 1:K)



p_guide_degas_xl <- pheatmap(jaccard_guide_degas_xl,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(0, 0.5, length.out = 101),
                             color = colorRampPalette(c("white", "firebrick3"))(100),
                             main = "Variant Jaccard Indices",
                             display_numbers = ifelse(jaccard_guide_degas_xl == 0, 0, round(jaccard_guide_degas_xl, 2)),
                             legend = FALSE,
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)

p_guide_degas_lt <- pheatmap(jaccard_guide_degas_lt,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(0, 0.5, length.out = 101),
                             color = colorRampPalette(c("white", "firebrick3"))(100),
                             main = "Trait Jaccard Indices",
                             display_numbers = ifelse(jaccard_guide_degas_lt == 0, 0, round(jaccard_guide_degas_lt, 2)),
                             legend = FALSE,
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)


## save at 19.5x11

grid.newpage()
grid.draw(
  arrangeGrob(
    nullGrob(),nullGrob(),nullGrob(), ## empty row
    p_guide_degas_xl$gtable, nullGrob(), p_guide_degas_lt$gtable,
    ncol = 3,
    nrow = 2,
    widths = c(1, 0.07, 1),
    heights = c(0.08, 1),   
    top = textGrob(
      "GUIDE & DeGAs Component Overlap",
      gp = gpar(fontsize = 35, fontface = "bold"),
      vjust = 1, hjust = 0.41           
    )
  )
)




p_guide_xl <- pheatmap(jaccard_guide_xl,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       breaks = seq(0, 0.8, length.out = 101),
                       color = colorRampPalette(c("white", "firebrick3"))(100),
                       main = "Variant Jaccard Indices",
                       display_numbers = ifelse(jaccard_guide_xl == 0, 0, round(jaccard_guide_xl, 2)),
                       legend = FALSE,
                       fontsize_number = 15,
                       fontsize_row = 16,
                       fontsize_col = 16,
                       fontsize = 20,
                       angle_col = 45)

p_guide_lt <- pheatmap(jaccard_guide_lt,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       breaks = seq(0, 0.8, length.out = 101),
                       color = colorRampPalette(c("white", "firebrick3"))(100),
                       main = "Trait Jaccard Indices",
                       display_numbers = ifelse(jaccard_guide_lt == 0, 0, round(jaccard_guide_lt, 2)),
                       legend = FALSE,
                       fontsize_number = 15,
                       fontsize_row = 16,
                       fontsize_col = 16,
                       fontsize = 20,
                       angle_col = 45)


## save at 19.5x11

grid.newpage()
grid.draw(
  arrangeGrob(
    nullGrob(),nullGrob(),nullGrob(), ## empty row
    p_guide_xl$gtable, nullGrob(), p_guide_lt$gtable,
    ncol = 3,
    nrow = 2,
    widths = c(1, 0.07, 1),
    heights = c(0.08, 1),   
    top = textGrob(
      "GUIDE Component Overlap",
      gp = gpar(fontsize = 35, fontface = "bold"),
      vjust = 1, hjust = 0.39            
    )
  )
)




p_degas_xl <- pheatmap(jaccard_degas_xl,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       breaks = seq(0, 0.8, length.out = 101),
                       color = colorRampPalette(c("white", "firebrick3"))(100),
                       main = "Variant Jaccard Indices",
                       display_numbers = ifelse(jaccard_degas_xl == 0, 0, round(jaccard_degas_xl, 2)),
                       legend = FALSE,
                       fontsize_number = 15,
                       fontsize_row = 16,
                       fontsize_col = 16,
                       fontsize = 20,
                       angle_col = 45)

p_degas_lt <- pheatmap(jaccard_degas_lt,
                       cluster_rows = FALSE,
                       cluster_cols = FALSE,
                       breaks = seq(0, 0.8, length.out = 101),
                       color = colorRampPalette(c("white", "firebrick3"))(100),
                       main = "Trait Jaccard Indices",
                       display_numbers = ifelse(jaccard_degas_lt == 0, 0, round(jaccard_degas_lt, 2)),
                       legend = FALSE,
                       fontsize_number = 15,
                       fontsize_row = 16,
                       fontsize_col = 16,
                       fontsize = 20,
                       angle_col = 45)


## save at 19.5x11

grid.newpage()
grid.draw(
  arrangeGrob(
    nullGrob(),nullGrob(),nullGrob(), ## empty row
    p_degas_xl$gtable, nullGrob(), p_degas_lt$gtable,
    ncol = 3,
    nrow = 2,
    widths = c(1, 0.07, 1),
    heights = c(0.08, 1),   
    top = textGrob(
      "DeGAs Component Overlap",
      gp = gpar(fontsize = 35, fontface = "bold"),
      vjust = 1, hjust = 0.38            
    )
  )
)





## Plots related to component sparsity


### Plotting cumulative contribution scores


plot_data_xl <- lapply(1:K, function(ii) {
  sorted_degas <- cumsum(sort(W.xl.degas[,ii]^2, decreasing = TRUE))
  sorted_guide <- cumsum(sort(W.xl.guide[,ii]^2, decreasing = TRUE))
  
  data.frame(
    Index = 0:length(sorted_guide),
    CumSum = c(0, sorted_guide),
    Method = "GUIDE",
    Group = paste0("GUIDE_", ii)
  ) %>%
    bind_rows(
      data.frame(
        Index = 0:length(sorted_degas),
        CumSum = c(0, sorted_degas),
        Method = "DeGAs",
        Group = paste0("DeGAs_", ii)
      )
    )
}) %>% bind_rows()

plot_xl = ggplot(plot_data_xl, aes(x = Index, y = CumSum, color = Method, group = Group)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Number of Variants",
    y = "Variant Cumulative Contribution Scores",
    title = "",
    fill = "Method"
  ) +
  scale_color_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +
  scale_y_continuous(breaks = seq(0,1,by=0.2)) +
  scale_x_continuous(breaks = seq(0,600,by=150)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(10, 80, 10, 80))





plot_data_lt <- lapply(1:K, function(ii) {
  sorted_degas <- cumsum(sort(W.lt.degas[,ii]^2, decreasing = TRUE))
  sorted_guide <- cumsum(sort(W.lt.guide[,ii]^2, decreasing = TRUE))
  
  data.frame(
    Index = 0:length(sorted_guide),
    CumSum = c(0, sorted_guide),
    Method = "GUIDE",
    Group = paste0("GUIDE_", ii)
  ) %>%
    bind_rows(
      data.frame(
        Index = 0:length(sorted_degas),
        CumSum = c(0, sorted_degas),
        Method = "DeGAs",
        Group = paste0("DeGAs_", ii)
      )
    )
}) %>% bind_rows()

plot_data_lt$Method <- factor(plot_data_lt$Method, levels = c("GUIDE", "DeGAs"))

plot_lt = ggplot(plot_data_lt, aes(x = Index, y = CumSum, color = Method, group = Group)) +
  geom_line(linewidth = 1.2) +
  labs(
    x = "Number of Traits",
    y = "Trait Cumulative Contribution Scores",
    title = "",
    fill = "Method"
  ) +
  scale_color_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +
  scale_y_continuous(breaks = seq(0,1,by=0.2)) +
  scale_x_continuous(breaks = seq(0,120,by=25)) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 23, vjust = -1),              
        axis.title.y = element_text(size = 23, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 10, 20, 15),
        legend.text = element_text(size = 16),  
        legend.title = element_text(size = 18),    
        legend.key.size = unit(1.7, "cm"))



combined_plot <- plot_grid(plot_xl, plot_lt, nrow = 1, align = "h")

title <- ggdraw() + 
  draw_label(
    "GUIDE & DeGAs Cumulative Contribution Scores", 
    size = 30, 
    hjust = 0.5,
    vjust = 1
  )

# save at 19.65x9
plot_grid(title, combined_plot, ncol = 1, rel_heights = c(0.05, 1))



## Plots for the kurtosis values for each component and the
## number of variants and traits loaded onto each, separately for GUIDE and DeGAs




k.guide.total = c()
k.degas.total = c()

for (ii in 1:K) {
  
  k.guide = round(kurtosis(c(W.xl.guide[,ii], W.lt.guide[,ii])) - 3, 1)
  k.degas = round(kurtosis(c(W.xl.degas[,ii], W.lt.degas[,ii])) - 3, 1)
  
  k.guide.total = c(k.guide.total, k.guide)
  k.degas.total = c(k.degas.total, k.degas)
}

cbind(k.guide.total, k.degas.total)


colSums(W.xl.guide_signif != 0)
colSums(W.xl.degas_signif != 0)

colSums(W.lt.guide_signif != 0)
colSums(W.lt.degas_signif != 0)


nrow(df) - sum(apply(W.xl.guide_signif, 1, function(row) all(row == 0)))
nrow(df) - sum(apply(W.xl.degas_signif, 1, function(row) all(row == 0)))

ncol(df) - sum(apply(W.lt.guide_signif, 1, function(row) all(row == 0)))
ncol(df) - sum(apply(W.lt.degas_signif, 1, function(row) all(row == 0)))




kurt.df <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), c(length(k.guide.total), length(k.degas.total))),
                  levels = c("GUIDE", "DeGAs")),
  Value = c(k.guide.total, k.degas.total)
)

kurt_plot = ggplot(kurt.df, aes(x = Method, y = Value, colour = Method)) +
  geom_beeswarm(cex = 1.8, size = 3) +
  scale_color_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0,350,by=50)) +
  labs(title = "Kurtosis Values",
       x = NULL,
       y = "Kurtosis" ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 3),
    axis.title.x = element_text(size = 19, vjust = -0.5),              
    axis.title.y = element_text(size = 19, vjust = 3),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),
    plot.margin = margin(20, 10, 20, 15),
    legend.position = "none"
  )





vars.load.df <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), c(K, K)),
                  levels = c("GUIDE", "DeGAs")),
  Value = c(colSums(W.xl.guide_signif != 0), colSums(W.xl.degas_signif != 0))
)

var_plot = ggplot(vars.load.df, aes(x = Method, y = Value, colour = Method)) +
  geom_beeswarm(cex = 1.3, size = 3) +
  scale_color_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0,180,by=20)) +
  labs(title = "Variants Loaded",
       x = NULL,
       y = "Number of Variants" ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 3),
    axis.title.x = element_text(size = 19, vjust = -0.5),              
    axis.title.y = element_text(size = 19, vjust = 3),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),
    plot.margin = margin(20, 10, 20, 15),
    legend.position = "none"
  )



traits.load.df <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), c(K, K)),
                  levels = c("GUIDE", "DeGAs")),
  Value = c(colSums(W.lt.guide_signif != 0), colSums(W.lt.degas_signif != 0))
)

trait_plot = ggplot(traits.load.df, aes(x = Method, y = Value, colour = Method)) +
  geom_beeswarm(cex = 1.3, size = 3) +
  scale_color_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0,40,by=10)) +
  labs(title = "Traits Loaded",
       x = NULL,
       y = "Number of Traits" ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5, vjust = 3),
    axis.title.x = element_text(size = 19, vjust = -0.5),              
    axis.title.y = element_text(size = 19, vjust = 3),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 16),
    plot.margin = margin(20, 10, 20, 15),
    legend.position = "right",
    legend.text = element_text(size = 16),       
    legend.title = element_text(size = 18, hjust = 0.5),   
    legend.key.size = unit(1.7, "cm")
  )


tight_theme <- theme(plot.margin = margin(0, -30, 0, -50))
plots_icasso_ica <- plot_grid(
  kurt_plot + tight_theme, var_plot + tight_theme, trait_plot + tight_theme, 
  ncol = 3, align = "hv"
)

icasso_ica_total <- ggdraw() +
  draw_label("GUIDE & DeGAs Component Sparsity Measures", fontface = 'plain', size = 28, hjust = 0.54, vjust = -20.5) +
  draw_plot(
    plot_grid(
      plots_icasso_ica + theme(plot.margin = margin(30, 30, 85, 35))
    ),
    y = -0.08,
    x = 0.025
  )

## save at 15.5x12.5 inches
icasso_ica_total









## Heatmaps for GUIDE and DeGAs variant and phenotype matrices (not included in thesis)



guide_heatmap_xl <- pheatmap(W.xl.guide,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(-1, 1, length.out = 101),
                             color = colorRampPalette(c("navy", "gainsboro", "firebrick3"))(100),
                             main = "Variant Weights",
                             legend = FALSE,
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)

guide_heatmap_lt <- pheatmap(W.lt.guide,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(-1, 1, length.out = 101),
                             color = colorRampPalette(c("navy", "gainsboro", "firebrick3"))(100),
                             main = "Trait Weights",
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)


## save at 17x9.5

grid.newpage()
grid.draw(
  arrangeGrob(
    nullGrob(),nullGrob(),nullGrob(), ## empty row
    guide_heatmap_xl$gtable, nullGrob(), guide_heatmap_lt$gtable,
    ncol = 3,
    nrow = 2,
    widths = c(1, 0.07, 1),
    heights = c(0.08, 1),   
    top = textGrob(
      "GUIDE Weights Heatmaps",
      gp = gpar(fontsize = 35, fontface = "bold"),
      vjust = 1, hjust = 0.5            
    )
  )
)






degas_heatmap_xl <- pheatmap(W.xl.degas,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(-1, 1, length.out = 101),
                             color = colorRampPalette(c("navy", "gainsboro", "firebrick3"))(100),
                             main = "Variant Weights",
                             legend = FALSE,
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)

degas_heatmap_lt <- pheatmap(W.lt.degas,
                             cluster_rows = FALSE,
                             cluster_cols = FALSE,
                             breaks = seq(-1, 1, length.out = 101),
                             color = colorRampPalette(c("navy", "gainsboro", "firebrick3"))(100),
                             main = "Trait Weights",
                             fontsize_number = 15,
                             fontsize_row = 16,
                             fontsize_col = 16,
                             fontsize = 20,
                             angle_col = 45)


## save at 19.5x11

grid.newpage()
grid.draw(
  arrangeGrob(
    nullGrob(),nullGrob(),nullGrob(), ## empty row
    degas_heatmap_xl$gtable, nullGrob(), degas_heatmap_lt$gtable,
    ncol = 3,
    nrow = 2,
    widths = c(1, 0.07, 1),
    heights = c(0.08, 1),   
    top = textGrob(
      "DeGAs Weights Heatmaps",
      gp = gpar(fontsize = 35, fontface = "bold"),
      vjust = 1, hjust = 0.5          
    )
  )
)

