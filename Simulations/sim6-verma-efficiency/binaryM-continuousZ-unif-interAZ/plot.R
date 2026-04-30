# packages
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)
library(stats)
library(xtable)
library(gridExtra)
library(cowplot)
library(here)

setwd(here::here("sim6-verma-efficiency/binaryM-continuousZ-unif-interAZ"))

# plot function
source("plot-sub.R")

# sample size
n.vec <- c(500,1000,2000,4000,8000)

# number of simulation
nsim <- 1000

### binary case ====
# the truth
load("../DGPs/6-truth-binaryM-continuousZ-unif-interAZ.Rdata")


## onestep====
load("Onestep/result.Rdata")
p.all.z.one <- plot.tmle(r'($\psi^{+}(\hat{Q};Z)$)',z="all.z")

# 
# p.tmle <- plot_grid(
#   p.z0.tmle,p.z1.tmle,p.all.z.tmle,p.opt.tmle,
#   align = "v", ncol = 1
# )
# 
# p.one <- plot_grid(
#   p.z0.one,p.z1.one,p.all.z.one,p.opt.one,
#   align = "v", ncol = 1
# )
# 
# 
# p.final <- plot_grid(
# p.one,
#   align = "h",
#   ncol = 1
# )


ggsave("plot.pdf", plot = p.all.z.one, width = 8, height = 6, units = "in")

