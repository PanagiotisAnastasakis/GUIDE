
library(ggplot2)
library(patchwork)
library(latex2exp)
library(cowplot)



results_misspec_recovery_guide = readRDS("results_misspec_recovery_guide.rds")

k.est = 20

get_quantiles <- function(values) {
  c(median(values), quantile(values, 0.025), quantile(values, 0.975))
}


## Gaussian sampler, low sparsity

## all components

quantiles.gaussian.low.misspec.rec = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["gaussian"]][["low"]][["all"]][,jj])))

quantiles.gaussian.low.misspec.rec = rbind(c(0,0,0), quantiles.gaussian.low.misspec.rec)

data.misspec.rec.gaussian.low = data.frame(
  k = 0:k.est,
  median = quantiles.gaussian.low.misspec.rec[,1],
  lower = quantiles.gaussian.low.misspec.rec[,2],
  upper = quantiles.gaussian.low.misspec.rec[,3]
)

## top 5 components

quantiles.gaussian.low.misspec.rec.top5 = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["gaussian"]][["low"]][["top5"]][,jj])))

quantiles.gaussian.low.misspec.rec.top5 = rbind(c(0,0,0), quantiles.gaussian.low.misspec.rec.top5)


data.misspec.rec.gaussian.low.top5 = data.frame(
  k = 0:k.est,
  median = quantiles.gaussian.low.misspec.rec.top5[,1],
  lower = quantiles.gaussian.low.misspec.rec.top5[,2],
  upper = quantiles.gaussian.low.misspec.rec.top5[,3]
)


## Gaussian sampler, high sparsity

## all components

quantiles.gaussian.high.misspec.rec = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["gaussian"]][["high"]][["all"]][,jj])))

quantiles.gaussian.high.misspec.rec = rbind(c(0,0,0), quantiles.gaussian.high.misspec.rec)

data.misspec.rec.gaussian.high = data.frame(
  k = 0:k.est,
  median = quantiles.gaussian.high.misspec.rec[,1],
  lower = quantiles.gaussian.high.misspec.rec[,2],
  upper = quantiles.gaussian.high.misspec.rec[,3]
)

## top 5 components

quantiles.gaussian.high.misspec.rec.top5 = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["gaussian"]][["high"]][["top5"]][,jj])))

quantiles.gaussian.high.misspec.rec.top5 = rbind(c(0,0,0), quantiles.gaussian.high.misspec.rec.top5)

data.misspec.rec.gaussian.high.top5 = data.frame(
  k = 0:k.est,
  median = quantiles.gaussian.high.misspec.rec.top5[,1],
  lower = quantiles.gaussian.high.misspec.rec.top5[,2],
  upper = quantiles.gaussian.high.misspec.rec.top5[,3]
)


## Laplace sampler, low sparsity

## all components

quantiles.laplace.low.misspec.rec = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["laplace"]][["low"]][["all"]][,jj])))

quantiles.laplace.low.misspec.rec = rbind(c(0,0,0), quantiles.laplace.low.misspec.rec)

data.misspec.rec.laplace.low = data.frame(
  k = 0:k.est,
  median = quantiles.laplace.low.misspec.rec[,1],
  lower = quantiles.laplace.low.misspec.rec[,2],
  upper = quantiles.laplace.low.misspec.rec[,3]
)

## top 5 components

quantiles.laplace.low.misspec.rec.top5 = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["laplace"]][["low"]][["top5"]][,jj])))

quantiles.laplace.low.misspec.rec.top5 = rbind(c(0,0,0), quantiles.laplace.low.misspec.rec.top5)

data.misspec.rec.laplace.low.top5 = data.frame(
  k = 0:k.est,
  median = quantiles.laplace.low.misspec.rec.top5[,1],
  lower = quantiles.laplace.low.misspec.rec.top5[,2],
  upper = quantiles.laplace.low.misspec.rec.top5[,3]
)


## Laplace sampler, high sparsity

## all components

quantiles.laplace.high.misspec.rec = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["laplace"]][["high"]][["all"]][,jj])))

quantiles.laplace.high.misspec.rec = rbind(c(0,0,0), quantiles.laplace.high.misspec.rec)

data.misspec.rec.laplace.high = data.frame(
  k = 0:k.est,
  median = quantiles.laplace.high.misspec.rec[,1],
  lower = quantiles.laplace.high.misspec.rec[,2],
  upper = quantiles.laplace.high.misspec.rec[,3]
)

## top 5 components

quantiles.laplace.high.misspec.rec.top5 = t(sapply(1:k.est, function(jj) 
  get_quantiles(results_misspec_recovery_guide[["laplace"]][["high"]][["top5"]][,jj])))

quantiles.laplace.high.misspec.rec.top5 = rbind(c(0,0,0), quantiles.laplace.high.misspec.rec.top5)

data.misspec.rec.laplace.high.top5 = data.frame(
  k = 0:k.est,
  median = quantiles.laplace.high.misspec.rec.top5[,1],
  lower = quantiles.laplace.high.misspec.rec.top5[,2],
  upper = quantiles.laplace.high.misspec.rec.top5[,3]
)






plot_rec <- function(data, title, y_breaks, yend = 10, vjust_y = 1) {
  ggplot(data, aes(x = k, y = median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#FF4DB8") + 
    geom_line(linewidth = 1, color = "#FF4DB8") +  
    annotate("segment", x = 0, y = 0, xend = 20, yend = yend, color = "blue", linewidth = 1, linetype = "dashed") +
    labs(title = title, x = "CQI-Ordered GUIDE Components", y = "Matched Components Found") +
    theme_minimal() +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_continuous(breaks = seq(0, 20, by = 4)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 20, hjust = 0.5),
          axis.title.x = element_text(size = 17, vjust = 0.2),
          axis.title.y = element_text(size = 17, vjust = vjust_y),
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 14),
          plot.margin = margin(20, 10, 20, 10)) +
    coord_cartesian(xlim = range(data$k), ylim = c(0, yend))
}


## all components plot

p.misspec.rec.gaussian.low <- plot_rec(data.misspec.rec.gaussian.low, "Gaussian Setting, Low Sparsity", y_breaks = seq(0, 10, by = 2))  
p.misspec.rec.gaussian.high <- plot_rec(data.misspec.rec.gaussian.high, "Gaussian Setting, High Sparsity", y_breaks = seq(0, 10, by = 2))  
p.misspec.rec.laplace.low <- plot_rec(data.misspec.rec.laplace.low, "Laplace Setting, Low Sparsity", y_breaks = seq(0, 10, by = 2))  
p.misspec.rec.laplace.high <- plot_rec(data.misspec.rec.laplace.high, "Laplace Setting, High Sparsity", y_breaks = seq(0, 10, by = 2))  

plots_gridmisspec_rec <- plot_grid(
  p.misspec.rec.gaussian.low, p.misspec.rec.gaussian.high,
  p.misspec.rec.laplace.low, p.misspec.rec.laplace.high,
  ncol = 2, align = "hv"
)


final_plot_gridmisspec_rec <- ggdraw() +
  draw_label("Matched Components Found with CQI Ordering", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_plot(
    plot_grid(
      plots_gridmisspec_rec + theme(plot.margin = margin(10, 20, 25, 15))),
    y = -0.03
  )


## save at 16x12 inches
final_plot_gridmisspec_rec





p.misspec.rec.gaussian.low.top5 <- plot_rec(data.misspec.rec.gaussian.low.top5, "Gaussian Sampler, Low Sparsity", y_breaks = 0:5, yend = 5, vjust_y = 3)  
p.misspec.rec.gaussian.high.top5 <- plot_rec(data.misspec.rec.gaussian.high.top5, "Gaussian Sampler, High Sparsity", y_breaks = 0:5, yend = 5, vjust_y = 3)  
p.misspec.rec.laplace.low.top5 <- plot_rec(data.misspec.rec.laplace.low.top5, "Laplace Sampler, Low Sparsity", y_breaks = 0:5, yend = 5, vjust_y = 3)  
p.misspec.rec.laplace.high.top5 <- plot_rec(data.misspec.rec.laplace.high.top5, "Laplace Sampler, High Sparsity", y_breaks = 0:5, yend = 5, vjust_y = 3)  

plots_gridmisspec_rec_top5 <- plot_grid(
  p.misspec.rec.gaussian.low.top5, p.misspec.rec.gaussian.high.top5,
  p.misspec.rec.laplace.low.top5, p.misspec.rec.laplace.high.top5,
  ncol = 2, align = "hv"
)

final_plot_gridmisspec_rec_top5 <- ggdraw() +
  draw_label("Top 5 True Components against CQI-Ordered GUIDE Components", fontface = 'plain', size = 28, hjust = 0.5, vjust = -20.3) +
  draw_plot(
    plot_grid(
      plots_gridmisspec_rec_top5 + theme(plot.margin = margin(10, 20, 25, 15))),
    y = -0.03
  )

## save at 16x12 inches
final_plot_gridmisspec_rec_top5








