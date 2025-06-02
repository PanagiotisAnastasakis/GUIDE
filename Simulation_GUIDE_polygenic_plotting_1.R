
source("/Users/panos/Desktop/MSc Thesis/Code/helper.R")
library(ggplot2)
library(patchwork)
library(latex2exp)
library(cowplot)
library(ExtDist)



results = readRDS("polygenic_results.rds")



#### Storing the data 


#### gaussian sampler across all configurations

n.reps = 300

## C1

data.c1.gaussian.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - None
    results[["C1"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C1"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C1"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - None
    results[["C1"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C1"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c1.gaussian.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C1"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C1"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C1"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C1"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C1"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c1.gaussian.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["gaussian"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C1"]][["gaussian"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C1"]][["gaussian"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C1"]][["gaussian"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C1"]][["gaussian"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C1"]][["gaussian"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)



## C2

data.c2.gaussian.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C2"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C2"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C2"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C2"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C2"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c2.gaussian.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C2"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C2"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C2"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C2"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C2"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c2.gaussian.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["gaussian"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C2"]][["gaussian"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C2"]][["gaussian"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C2"]][["gaussian"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C2"]][["gaussian"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C2"]][["gaussian"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)



## C3

data.c3.gaussian.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C3"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C3"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C3"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C3"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C3"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c3.gaussian.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C3"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C3"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C3"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C3"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C3"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c3.gaussian.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["gaussian"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C3"]][["gaussian"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C3"]][["gaussian"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C3"]][["gaussian"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C3"]][["gaussian"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C3"]][["gaussian"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)




## C4

data.c4.gaussian.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C4"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C4"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C4"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C4"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C4"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c4.gaussian.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["gaussian"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C4"]][["gaussian"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C4"]][["gaussian"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C4"]][["gaussian"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C4"]][["gaussian"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C4"]][["gaussian"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c4.gaussian.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["gaussian"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C4"]][["gaussian"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C4"]][["gaussian"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C4"]][["gaussian"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C4"]][["gaussian"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C4"]][["gaussian"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)




#### laplace sampler across all configurations



## C1

data.c1.laplace.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C1"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C1"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C1"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C1"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C1"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c1.laplace.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C1"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C1"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C1"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C1"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C1"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c1.laplace.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C1"]][["laplace"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C1"]][["laplace"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C1"]][["laplace"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C1"]][["laplace"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C1"]][["laplace"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C1"]][["laplace"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)



## C2

data.c2.laplace.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C2"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C2"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C2"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C2"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C2"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c2.laplace.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C2"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C2"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C2"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C2"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C2"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c2.laplace.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C2"]][["laplace"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C2"]][["laplace"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C2"]][["laplace"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C2"]][["laplace"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C2"]][["laplace"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C2"]][["laplace"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)



## C3

data.c3.laplace.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C3"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C3"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C3"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C3"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C3"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c3.laplace.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C3"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C3"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C3"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C3"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C3"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c3.laplace.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C3"]][["laplace"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C3"]][["laplace"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C3"]][["laplace"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C3"]][["laplace"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C3"]][["laplace"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C3"]][["laplace"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)




## C4

data.c4.laplace.xl <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Low
    results[["C4"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - Medium
    results[["C4"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.xl,  # GUIDE - High
    results[["C4"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Low
    results[["C4"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.xl, # DeGAs - Medium
    results[["C4"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.xl  # DeGAs - High
  )
)


data.c4.laplace.lt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["laplace"]][["none"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Low
    results[["C4"]][["laplace"]][["low"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - Medium
    results[["C4"]][["laplace"]][["high"]][["GUIDE"]]$avg.cor.lt,  # GUIDE - High
    results[["C4"]][["laplace"]][["none"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Low
    results[["C4"]][["laplace"]][["low"]][["DeGAs"]]$avg.cor.lt, # DeGAs - Medium
    results[["C4"]][["laplace"]][["high"]][["DeGAs"]]$avg.cor.lt  # DeGAs - High
  )
)


data.c4.laplace.kurt <- data.frame(
  Method = factor(rep(c("GUIDE", "DeGAs"), each = 3*n.reps), levels = c("GUIDE", "DeGAs")),
  Setting = factor(rep(rep(c("None", "Low", "High"), each = n.reps), 2), levels = c("None", "Low", "High")),
  Value = c(
    results[["C4"]][["laplace"]][["none"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Low
    results[["C4"]][["laplace"]][["low"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - Medium
    results[["C4"]][["laplace"]][["high"]][["GUIDE"]]$avg.kurtosis,  # GUIDE - High
    results[["C4"]][["laplace"]][["none"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Low
    results[["C4"]][["laplace"]][["low"]][["DeGAs"]]$avg.kurtosis, # DeGAs - Medium
    results[["C4"]][["laplace"]][["high"]][["DeGAs"]]$avg.kurtosis  # DeGAs - High
  )
)





#### Plotting


## Figure 1, Gaussian setting for C1, C2, C3, C5

## C1

f.c1.gaussian.xl = ggplot(data.c1.gaussian.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) + 
  scale_y_continuous(breaks = c(0.4,0.5, 0.6, 0.7, 0.8, 0.9)) +
  theme(
    legend.position = "none",           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c1.gaussian.lt = ggplot(data.c1.gaussian.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(#title = "Simulation Results for Gaussian Sampler",
       title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5, vjust = 1.5),         
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
    )




f.c1.gaussian.kurt = ggplot(data.c1.gaussian.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )


## Save at 16 x 4.85 inches

#plot_grid(f.c1.gaussian.xl, f.c1.gaussian.lt, f.c1.gaussian.kurt, nrow = 1, align = "h", labels = c('C1'), label_size = 25)
plot_grid(f.c1.gaussian.xl, f.c1.gaussian.lt, f.c1.gaussian.kurt, nrow = 1, align = "h", labels = c('C1'), label_size = 25) +
  draw_label("Simulation Results with Gaussian Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -7.65)



## C2

f.c2.gaussian.xl = ggplot(data.c2.gaussian.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +  # Left-align title +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) +
  theme(
    plot.title = element_text(size = 28), ## we must have this line at least in one plot in each row to make "C" label be higher
    legend.position = "none",           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c2.gaussian.lt = ggplot(data.c2.gaussian.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "none",  
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )




f.c2.gaussian.kurt = ggplot(data.c2.gaussian.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",             
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )


#plot_grid(f.c2.gaussian.xl, f.c2.gaussian.lt, f.c2.gaussian.kurt, nrow = 1, align = "h", labels = c('C2'), label_size = 25)

plot_grid(f.c2.gaussian.xl, f.c2.gaussian.lt, f.c2.gaussian.kurt, nrow = 1, align = "h",
           labels = c('C1'), label_size = 25) +
  draw_label("Simulation Results with Gaussian Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -10)

## save at 16x5

ggdraw() +
  draw_plot(
    plot_grid(f.c2.gaussian.xl, f.c2.gaussian.lt, f.c2.gaussian.kurt,
              nrow = 1, align = "h", labels = c('C1'), label_size = 25),
    y = -0.038  # Push down by lowering y (0 = bottom, 1 = top)
  ) +
  draw_label("Simulation Results with Gaussian Setting",
             fontface = 'plain', size = 28, hjust = 0.5, y = 0.96)





## C3

f.c3.gaussian.xl = ggplot(data.c3.gaussian.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +  # Left-align title +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    plot.title = element_text(size = 28), ## we must have this line at least in one plot in each row to make "C" label be higher
    legend.position = "none",             
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



# f.c3.gaussian.lt = ggplot(data.c3.gaussian.lt, aes(x = Setting, y = Value, fill = Method)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
#   theme_minimal() +
#   labs(title = "",
#        x = "",
#        y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
#        fill = "Method") +
#   scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
#   theme(
#     legend.position = "none",        
#     axis.title.y = element_text(size = 17),
#     axis.text.x = element_text(size = 16),  
#     axis.text.y = element_text(size = 14)
#   )

f.c3.gaussian.lt = ggplot(data.c3.gaussian.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "bottom",          
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),        # Make legend text bigger
    legend.title = element_text(size = 18),       # Make legend title bigger
    legend.key.size = unit(1.7, "cm")  
  )






f.c3.gaussian.kurt = ggplot(data.c3.gaussian.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",         
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )

## save at 16x6.1
plot_grid(f.c3.gaussian.xl, f.c3.gaussian.lt, f.c3.gaussian.kurt, nrow = 1, align = "h", labels = c('C2'), label_size = 25)


## C4

f.c4.gaussian.xl = ggplot(data.c4.gaussian.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
    x = "",
    y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
    fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28),   ## we must have this line at least in one plot in each row to make "C" label be higher
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c4.gaussian.lt = ggplot(data.c4.gaussian.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "bottom",          
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),        # Make legend text bigger
    legend.title = element_text(size = 18),       # Make legend title bigger
    legend.key.size = unit(1.7, "cm")  
  )




f.c4.gaussian.kurt = ggplot(data.c4.gaussian.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",             
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)  
  )


## Save at 16x5.1 inches

plot_grid(f.c4.gaussian.xl, f.c4.gaussian.lt, f.c4.gaussian.kurt, nrow = 1, align = "h", labels = c('C4'), label_size = 25)






## Figure 2, Laplace sampler for C1, C2, C3, C5

## C1

f.c1.laplace.xl = ggplot(data.c1.laplace.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +  # Left-align title +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)) +
  theme(
    legend.position = "none",        
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c1.laplace.lt = ggplot(data.c1.laplace.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28, hjust = 0.5, vjust = 1.5),           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)         # Make legend key (box) bigger
    
  )




f.c1.laplace.kurt = ggplot(data.c1.laplace.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",            
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )


## Save at 16 x 4.85 inches

plot_grid(f.c1.laplace.xl, f.c1.laplace.lt, f.c1.laplace.kurt, nrow = 1, align = "h", labels = c('C1'), label_size = 25) +
  draw_label("Simulation Results with Laplace Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -7.65)



## C2

f.c2.laplace.xl = ggplot(data.c2.laplace.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +  # Left-align title +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28),  ## we must have this line at one plot in each row to make "C" label be higher          
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c2.laplace.lt = ggplot(data.c2.laplace.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )




f.c2.laplace.kurt = ggplot(data.c2.laplace.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",         
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )



#plot_grid(f.c2.laplace.xl, f.c2.laplace.lt, f.c2.laplace.kurt, nrow = 1, align = "h", labels = c('C2'), label_size = 25)

plot_grid(f.c2.laplace.xl, f.c2.laplace.lt, f.c2.laplace.kurt, nrow = 1, align = "h", labels = c('C1'), label_size = 25) +
  draw_label("Simulation Results with Laplace Setting", fontface = 'plain', size = 28, hjust = 0.5, vjust = -7.65)


ggdraw() +
  draw_plot(
    plot_grid(f.c2.laplace.xl, f.c2.laplace.lt, f.c2.laplace.kurt,
              nrow = 1, align = "h", labels = c('C1'), label_size = 25),
    y = -0.038  # Push down by lowering y (0 = bottom, 1 = top)
  ) +
  draw_label("Simulation Results with Laplace Setting",
             fontface = 'plain', size = 28, hjust = 0.5, y = 0.96)






## C3

f.c3.laplace.xl = ggplot(data.c3.laplace.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") +  # Left-align title +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28),  ## we must have this line in one plot in each row to make "C" label be higher               
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



# f.c3.laplace.lt = ggplot(data.c3.laplace.lt, aes(x = Setting, y = Value, fill = Method)) +
#   geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
#   theme_minimal() +
#   labs(title = "",
#        x = "",
#        y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
#        fill = "Method") +
#   scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
#   theme(
#     legend.position = "none",          
#     axis.title.y = element_text(size = 17),
#     axis.text.x = element_text(size = 16),  
#     axis.text.y = element_text(size = 14)
#   )


f.c3.laplace.lt = ggplot(data.c3.laplace.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "bottom",       
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),        # Make legend text bigger
    legend.title = element_text(size = 18),       # Make legend title bigger
    legend.key.size = unit(1.7, "cm")  
  )




f.c3.laplace.kurt = ggplot(data.c3.laplace.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",           
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)
  )


plot_grid(f.c3.laplace.xl, f.c3.laplace.lt, f.c3.laplace.kurt, nrow = 1, align = "h", labels = c('C2'), label_size = 25)


## C4

f.c4.laplace.xl = ggplot(data.c4.laplace.xl, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{XL}, \tilde{\textbf{W}}_{XL}$))"),
       fill = "Method") + 
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",
    plot.title = element_text(size = 28),  ## we must have this line in one plot in each row to make "C" label be higher              
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14) 
  )



f.c4.laplace.lt = ggplot(data.c4.laplace.lt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(Cor($\textbf{W}_{TL}, \tilde{\textbf{W}}_{TL}$))"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  scale_y_continuous(breaks = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
  theme(
    legend.position = "bottom",       
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 16),        # Make legend text bigger
    legend.title = element_text(size = 18),       # Make legend title bigger
    legend.key.size = unit(1.7, "cm")  
  )




f.c4.laplace.kurt = ggplot(data.c4.laplace.kurt, aes(x = Setting, y = Value, fill = Method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Dodge to separate groups
  theme_minimal() +
  labs(title = "",
       x = "",
       y = TeX(r"(\kappa$(\tilde{\textbf{W}}_{XL} \, | \, \tilde{\textbf{W}}_{TL})$)"),
       fill = "Method") +
  scale_fill_manual(values = c("GUIDE" = "#E69F00", "DeGAs" = "#56B4E9")) +  # Custom Colors
  theme(
    legend.position = "none",            
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16),  
    axis.text.y = element_text(size = 14)  
  )


## Save at 16x5.1 inches

plot_grid(f.c4.laplace.xl, f.c4.laplace.lt, f.c4.laplace.kurt, nrow = 1, align = "h", labels = c('C4'), label_size = 25)






### Investigating what is happening with the outliers

set.seed(4)

sampler = "laplace"
sparsity = "low"
n = 3000
p = 150
k = 20


sim_polygenic = polygenic_sim(n = n, p = p, k = k, sampler = sampler, sparsity = sparsity, return.afs = T)

W = sim_polygenic$W
W.xl = sim_polygenic$W.xl
W.lt = sim_polygenic$W.lt
afs = sim_polygenic$afs

weights.true = rbind(W.xl, W.lt)

ica_clustering = get_ica_clustering(W = W, K = k, reps = 20)
a = get_optimal_unmixing_matrix(ica_clustering$unmix.total, ica_clustering$clusters, ica_clustering$cors.unmix)

guide.list = get_guide(W, K = k, unmixing.matrix = a, verbose = F)

W.xl.guide = guide.list$W.xl
W.lt.guide = guide.list$W.lt
weights.guide = rbind(W.xl.guide, W.lt.guide)

weights.guide.matched = get_matched_weights(weights.true = weights.true,
                                            weights.estimated = weights.guide,
                                            n.variants = n)

W.xl.guide.matched = weights.guide.matched$weights.matched.xl
W.lt.guide.matched = weights.guide.matched$weights.matched.lt


#sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.guide.matched[,ii]))
#sapply(1:k, function(ii) cor(W.lt[,ii], W.lt.guide.matched[,ii]))




## save at 12.5, 9
plot_simulation(W.xl, W.xl.guide.matched, xlab = "GUIDE Variant Weights", ylab = "Simulated Variant Weights", main = "GUIDE & Simulated Variant Weights")
plot_simulation(W.lt, W.lt.guide.matched, xlab = "GUIDE Trait Weights", ylab = "Simulated Trait Weights", main = "GUIDE & Simulated Trait Weights")

#sapply(1:k, function(ii) cor(W.xl[top_afs_idx,][,ii], W.xl.guide.matched[top_afs_idx,][,ii]))
#sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.guide.matched[,ii]))

bottom_n_percent_indices <- function(n, vec) {
  
  threshold <- quantile(vec, n)
  which(vec <= threshold)
}

ind = 7
perc = 0.05


W.xl.guide.comp = W.xl.guide.matched[,ind]
W.xl.comp = W.xl[,ind]
top_afs_idx = bottom_n_percent_indices(perc, afs)



df <- data.frame(W.m = W.xl.guide.comp, W = W.xl.comp)

df_subset <- df[top_afs_idx, ]  

ggplot() +
  geom_point(data = df, aes(x = W.m, y = W), color = "black", shape = 16) +
  geom_point(data = df_subset, aes(x = W.m, y = W), color = "red", shape = 4, size = 4) +
  labs(x = "GUIDE Variant Weights", y = "True Variant Weights", title = "GUIDE & Simulated Variant Weights for one Component") +
  theme_minimal() +
  theme(plot.title = element_text(size = 28, hjust = 0.5, vjust = 3),
        axis.title.x = element_text(size = 22, vjust = -3),
        axis.title.y = element_text(size = 22, vjust = 3),
        axis.text.x = element_text(size = 16),  
        axis.text.y = element_text(size = 16),
        plot.margin = margin(20, 15, 25, 15))


#sapply(1:k, function(ii) cor(W.xl[top_afs_idx,][,ii], W.xl.guide.matched[top_afs_idx,][,ii]))
#sapply(1:k, function(ii) cor(W.xl[,ii], W.xl.guide.matched[,ii]))

