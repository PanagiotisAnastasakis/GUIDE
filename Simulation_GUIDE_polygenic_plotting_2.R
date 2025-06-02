
library(ggplot2)
library(patchwork)
library(latex2exp)
library(cowplot)

results_misspec = readRDS("polygenic_results_misspec.rds")


get_quantiles <- function(values) {
  c(median(values), quantile(values, 0.025), quantile(values, 0.975))
}

k.est.vals = 2:30

## Gaussian, low sparsity

quantiles.gaussian.low.guide.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.xl)))

quantiles.gaussian.low.guide.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.lt)))

quantiles.gaussian.low.guide.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.total)))

quantiles.gaussian.low.degas.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.xl)))

quantiles.gaussian.low.degas.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.lt)))

quantiles.gaussian.low.degas.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.total)))


data.misspec.gaussian.low.xl = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.low.guide.xl[,1], quantiles.gaussian.low.degas.xl[,1]),
  lower = c(quantiles.gaussian.low.guide.xl[,2], quantiles.gaussian.low.degas.xl[,2]),
  upper = c(quantiles.gaussian.low.guide.xl[,3], quantiles.gaussian.low.degas.xl[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.gaussian.low.lt = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.low.guide.lt[,1], quantiles.gaussian.low.degas.lt[,1]),
  lower = c(quantiles.gaussian.low.guide.lt[,2], quantiles.gaussian.low.degas.lt[,2]),
  upper = c(quantiles.gaussian.low.guide.lt[,3], quantiles.gaussian.low.degas.lt[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.gaussian.low.total = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.low.guide.total[,1], quantiles.gaussian.low.degas.total[,1]),
  lower = c(quantiles.gaussian.low.guide.total[,2], quantiles.gaussian.low.degas.total[,2]),
  upper = c(quantiles.gaussian.low.guide.total[,3], quantiles.gaussian.low.degas.total[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)


quantiles.gaussian.low.cor.min.guide = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["low"]][[k.est]][["GUIDE"]]$cor.min)))


data.misspec.gaussian.low.cor.min.guide = data.frame(
  k = k.est.vals,
  median = quantiles.gaussian.low.cor.min.guide[,1],
  lower = quantiles.gaussian.low.cor.min.guide[,2],
  upper = quantiles.gaussian.low.cor.min.guide[,3]
)



## Gaussian, high sparsity

quantiles.gaussian.high.guide.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.xl)))

quantiles.gaussian.high.guide.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.lt)))

quantiles.gaussian.high.guide.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.total)))

quantiles.gaussian.high.degas.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.xl)))

quantiles.gaussian.high.degas.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.lt)))

quantiles.gaussian.high.degas.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.total)))


data.misspec.gaussian.high.xl = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.high.guide.xl[,1], quantiles.gaussian.high.degas.xl[,1]),
  lower = c(quantiles.gaussian.high.guide.xl[,2], quantiles.gaussian.high.degas.xl[,2]),
  upper = c(quantiles.gaussian.high.guide.xl[,3], quantiles.gaussian.high.degas.xl[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.gaussian.high.lt = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.high.guide.lt[,1], quantiles.gaussian.high.degas.lt[,1]),
  lower = c(quantiles.gaussian.high.guide.lt[,2], quantiles.gaussian.high.degas.lt[,2]),
  upper = c(quantiles.gaussian.high.guide.lt[,3], quantiles.gaussian.high.degas.lt[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.gaussian.high.total = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.gaussian.high.guide.total[,1], quantiles.gaussian.high.degas.total[,1]),
  lower = c(quantiles.gaussian.high.guide.total[,2], quantiles.gaussian.high.degas.total[,2]),
  upper = c(quantiles.gaussian.high.guide.total[,3], quantiles.gaussian.high.degas.total[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)


quantiles.gaussian.high.cor.min.guide = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["gaussian"]][["high"]][[k.est]][["GUIDE"]]$cor.min)))


data.misspec.gaussian.high.cor.min.guide = data.frame(
  k = k.est.vals,
  median = quantiles.gaussian.high.cor.min.guide[,1],
  lower = quantiles.gaussian.high.cor.min.guide[,2],
  upper = quantiles.gaussian.high.cor.min.guide[,3]
)



## Laplace, low sparsity

quantiles.laplace.low.guide.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.xl)))

quantiles.laplace.low.guide.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.lt)))

quantiles.laplace.low.guide.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["GUIDE"]]$avg.cor.total)))

quantiles.laplace.low.degas.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.xl)))

quantiles.laplace.low.degas.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.lt)))

quantiles.laplace.low.degas.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["DeGAs"]]$avg.cor.total)))


data.misspec.laplace.low.xl = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.low.guide.xl[,1], quantiles.laplace.low.degas.xl[,1]),
  lower = c(quantiles.laplace.low.guide.xl[,2], quantiles.laplace.low.degas.xl[,2]),
  upper = c(quantiles.laplace.low.guide.xl[,3], quantiles.laplace.low.degas.xl[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.laplace.low.lt = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.low.guide.lt[,1], quantiles.laplace.low.degas.lt[,1]),
  lower = c(quantiles.laplace.low.guide.lt[,2], quantiles.laplace.low.degas.lt[,2]),
  upper = c(quantiles.laplace.low.guide.lt[,3], quantiles.laplace.low.degas.lt[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.laplace.low.total = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.low.guide.total[,1], quantiles.laplace.low.degas.total[,1]),
  lower = c(quantiles.laplace.low.guide.total[,2], quantiles.laplace.low.degas.total[,2]),
  upper = c(quantiles.laplace.low.guide.total[,3], quantiles.laplace.low.degas.total[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)


quantiles.laplace.low.cor.min.guide = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["low"]][[k.est]][["GUIDE"]]$cor.min)))


data.misspec.laplace.low.cor.min.guide = data.frame(
  k = k.est.vals,
  median = quantiles.laplace.low.cor.min.guide[,1],
  lower = quantiles.laplace.low.cor.min.guide[,2],
  upper = quantiles.laplace.low.cor.min.guide[,3]
)



## Laplace, high sparsity

quantiles.laplace.high.guide.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.xl)))

quantiles.laplace.high.guide.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.lt)))

quantiles.laplace.high.guide.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["GUIDE"]]$avg.cor.total)))

quantiles.laplace.high.degas.xl = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.xl)))

quantiles.laplace.high.degas.lt = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.lt)))

quantiles.laplace.high.degas.total = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["DeGAs"]]$avg.cor.total)))


data.misspec.laplace.high.xl = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.high.guide.xl[,1], quantiles.laplace.high.degas.xl[,1]),
  lower = c(quantiles.laplace.high.guide.xl[,2], quantiles.laplace.high.degas.xl[,2]),
  upper = c(quantiles.laplace.high.guide.xl[,3], quantiles.laplace.high.degas.xl[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.laplace.high.lt = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.high.guide.lt[,1], quantiles.laplace.high.degas.lt[,1]),
  lower = c(quantiles.laplace.high.guide.lt[,2], quantiles.laplace.high.degas.lt[,2]),
  upper = c(quantiles.laplace.high.guide.lt[,3], quantiles.laplace.high.degas.lt[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)

data.misspec.laplace.high.total = data.frame(
  k = rep(k.est.vals, 2),
  median = c(quantiles.laplace.high.guide.total[,1], quantiles.laplace.high.degas.total[,1]),
  lower = c(quantiles.laplace.high.guide.total[,2], quantiles.laplace.high.degas.total[,2]),
  upper = c(quantiles.laplace.high.guide.total[,3], quantiles.laplace.high.degas.total[,3]),
  Method = factor(rep(c("GUIDE", "DeGAs"), each = length(k.est.vals)), levels = c("GUIDE", "DeGAs"))
)


quantiles.laplace.high.cor.min.guide = t(sapply(k.est.vals, function(k.est) 
  get_quantiles(results_misspec[["laplace"]][["high"]][[k.est]][["GUIDE"]]$cor.min)))


data.misspec.laplace.high.cor.min.guide = data.frame(
  k = k.est.vals,
  median = quantiles.laplace.high.cor.min.guide[,1],
  lower = quantiles.laplace.high.cor.min.guide[,2],
  upper = quantiles.laplace.high.cor.min.guide[,3]
)


### Plotting


colors <- c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")

plot_misspec <- function(data, title, ylabel) {
  ggplot(data, aes(x = k, y = median, color = Method, fill = Method)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(title = title, x = "Number of Components", y = ylabel) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 17),
          axis.title.y = element_text(size = 17, vjust = 3),
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 14)) +
    xlim(0, 30) + 
    ylim(0.3, 1)
}

## Plotting X->L correlations

p.misspec.gaussian.low.xl <- plot_misspec(data.misspec.gaussian.low.xl, "", TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"))
p.misspec.gaussian.low.lt <- plot_misspec(data.misspec.gaussian.low.lt, "", TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"))  
p.misspec.gaussian.high.xl <- plot_misspec(data.misspec.gaussian.high.xl, "", TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"))  
p.misspec.gaussian.high.lt <- plot_misspec(data.misspec.gaussian.high.lt, "", TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"))  

legend <- suppressWarnings(get_legend(p.misspec.gaussian.low.xl + theme(legend.position = "right",
                                                                        legend.text = element_text(size = 16),        # Make legend text bigger
                                                                        legend.title = element_text(size = 18),       # Make legend title bigger
                                                                        legend.key.size = unit(1.7, "cm"))))

p.misspec.gaussian.low.xl <- p.misspec.gaussian.low.xl + theme(legend.position = "none")
p.misspec.gaussian.low.lt <- p.misspec.gaussian.low.lt + theme(legend.position = "none")
p.misspec.gaussian.high.xl <- p.misspec.gaussian.high.xl + theme(legend.position = "none")
p.misspec.gaussian.high.lt <- p.misspec.gaussian.high.lt + theme(legend.position = "none")

plots_grid_misspec_gaussian <- plot_grid(
  p.misspec.gaussian.low.xl, p.misspec.gaussian.low.lt,
  p.misspec.gaussian.high.xl, p.misspec.gaussian.high.lt,
  ncol = 2, align = "hv"
)


final_plot_misspec_gaussian <- ggdraw() +
  draw_label("Model Misspecification Simulation with Gaussian Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_label("Low Sparsity", x = 0.065, y = 0.96, size = 20, fontface = "bold") +
  draw_label("High Sparsity", x = 0.065, y = 0.48, size = 20, fontface = "bold") +
  draw_plot(
    plot_grid(
      plots_grid_misspec_gaussian + theme(plot.margin = margin(10, 10, 25, 15)), legend,
      ncol = 2, rel_widths = c(1, 0.11)
    ),
    y = -0.03
  )

## save at 16x12 inches
final_plot_misspec_gaussian


## Plotting Laplace

p.misspec.laplace.low.xl <- plot_misspec(data.misspec.laplace.low.xl, "", TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"))  
p.misspec.laplace.low.lt <- plot_misspec(data.misspec.laplace.low.lt, "", TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"))  
p.misspec.laplace.high.xl <- plot_misspec(data.misspec.laplace.high.xl, "", TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"))  
p.misspec.laplace.high.lt <- plot_misspec(data.misspec.laplace.high.lt, "", TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"))  

legend <- suppressWarnings(get_legend(p.misspec.laplace.low.xl + theme(legend.position = "right",
                                                                       legend.text = element_text(size = 16),        # Make legend text bigger
                                                                       legend.title = element_text(size = 18),       # Make legend title bigger
                                                                       legend.key.size = unit(1.7, "cm"))))

p.misspec.laplace.low.xl <- p.misspec.laplace.low.xl + theme(legend.position = "none")
p.misspec.laplace.low.lt <- p.misspec.laplace.low.lt + theme(legend.position = "none")
p.misspec.laplace.high.xl <- p.misspec.laplace.high.xl + theme(legend.position = "none")
p.misspec.laplace.high.lt <- p.misspec.laplace.high.lt + theme(legend.position = "none")

plots_grid_misspec_laplace <- plot_grid(
  p.misspec.laplace.low.xl, p.misspec.laplace.low.lt,
  p.misspec.laplace.high.xl, p.misspec.laplace.high.lt,
  ncol = 2, align = "hv"
)

final_plot_misspec_laplace <- ggdraw() +
  draw_label("Model Misspecification Simulation with Laplace Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_label("Low Sparsity", x = 0.065, y = 0.96, size = 20, fontface = "bold") +
  draw_label("High Sparsity", x = 0.065, y = 0.48, size = 20, fontface = "bold") +
  draw_plot(
    plot_grid(
      plots_grid_misspec_laplace + theme(plot.margin = margin(10, 10, 25, 15)), legend,
      ncol = 2, rel_widths = c(1, 0.12)
    ),
    y = -0.03
  )

## save at 16x12 inches
final_plot_misspec_laplace





##### Not used








## Plotting total correlations

p.misspec.gaussian.low.total <- plot_misspec(data.misspec.gaussian.low.total, "Gaussian Sampler, Low Sparsity", "Total Correlation")  
p.misspec.gaussian.high.total <- plot_misspec(data.misspec.gaussian.high.total, "Gaussian Sampler, High Sparsity", "Total Correlation")  
p.misspec.laplace.low.total <- plot_misspec(data.misspec.laplace.low.total, "Laplace Sampler, Low Sparsity", "Total Correlation")  
p.misspec.laplace.high.total <- plot_misspec(data.misspec.laplace.high.total, "Laplace Sampler, High Sparsity", "Total Correlation")  

legend <- suppressWarnings(get_legend(p.misspec.gaussian.low.total + theme(legend.position = "right",
                                                                       legend.text = element_text(size = 16),        # Make legend text bigger
                                                                       legend.title = element_text(size = 18),       # Make legend title bigger
                                                                       legend.key.size = unit(1.7, "cm"))))

p.misspec.gaussian.low.total <- p.misspec.gaussian.low.total + theme(legend.position = "none",
                                                                     plot.title = element_text(size = 20, hjust = 0.5))
p.misspec.gaussian.high.total <- p.misspec.gaussian.high.total + theme(legend.position = "none",
                                                                       plot.title = element_text(size = 20, hjust = 0.5))
p.misspec.laplace.low.total <- p.misspec.laplace.low.total + theme(legend.position = "none",
                                                                   plot.title = element_text(size = 20, hjust = 0.5))
p.misspec.laplace.high.total <- p.misspec.laplace.high.total + theme(legend.position = "none",
                                                                     plot.title = element_text(size = 20, hjust = 0.5))

plots_grid_misspec_total <- plot_grid(
  p.misspec.gaussian.low.total, p.misspec.gaussian.high.total,
  p.misspec.laplace.low.total, p.misspec.laplace.high.total,
  ncol = 2, align = "hv"
)

final_plot_misspec_total <- ggdraw() +
  draw_label("Model Misspecification Simulation Total Correlations", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_plot(
    plot_grid(
      plots_grid_misspec_total + theme(plot.margin = margin(10, 10, 25, 15)), legend,
      ncol = 2, rel_widths = c(1, 0.12)
    ),
    y = -0.03
  )

## save at 16x12 inches
final_plot_misspec_total










### plotting minimum correlation


plot_min <- function(data, title) {
  ggplot(data, aes(x = k, y = median, color = "#D55E00")) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, fill = "#D55E00") +
    geom_line(linewidth = 1, color = "#D55E00") +
    labs(title = title, x = "Number of Components",
         y = "Minimum Correlation") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title.x = element_text(size = 17, vjust = 0.2),
          axis.title.y = element_text(size = 17, vjust = 3),
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 14),
          plot.margin = margin(20, 10, 20, 10)) +
    xlim(0, 30) + 
    ylim(0, 1)
}


p.misspec.gaussian.low.cor.min.guide <- plot_min(data.misspec.gaussian.low.cor.min.guide, "Gaussian Sampler, Low Sparsity")  
p.misspec.gaussian.high.cor.min.guide <- plot_min(data.misspec.gaussian.high.cor.min.guide, "Gaussian Sampler, High Sparsity")  
p.misspec.laplace.low.cor.min.guide <- plot_min(data.misspec.laplace.low.cor.min.guide, "Laplace Sampler, Low Sparsity")  
p.misspec.laplace.high.cor.min.guide <- plot_min(data.misspec.laplace.high.cor.min.guide, "Laplace Sampler, High Sparsity")  

plots_grid_misspec_min <- plot_grid(
  p.misspec.gaussian.low.cor.min.guide, p.misspec.gaussian.high.cor.min.guide,
  p.misspec.laplace.low.cor.min.guide, p.misspec.laplace.high.cor.min.guide,
  ncol = 2, align = "hv"
)


final_plot_misspec_min <- ggdraw() +
  draw_label("Minimum Correlation of GUIDE Components under Model Misspecification", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_plot(
    plot_grid(
      plots_grid_misspec_min + theme(plot.margin = margin(10, 20, 25, 15))),
    y = -0.03
  )

## save at 16x12 inches
final_plot_misspec_min









