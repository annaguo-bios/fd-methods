library(here)
setwd(here("./sim1-consistency/multiM-d4/"))

# packages
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)
library(stats)
library(xtable)
library(gridExtra)
library(cowplot)

# plot function
source("plot-sub.R")

# sample size
n.vec <- c(250,500,1000,2000,4000,8000)

# number of simulation
nsim <- 1000

################################################ ATE ###################################################################
### multiM-d4 case ====
# the truth
load("../DGPs/1-truth-multi-d4.Rdata")

## 2D est2 ====
load("TMLE-est2a/result.Rdata")
p.2D.est2 <- plot.tmle(r'($\psi_{2a}(\hat{Q}^*)$)',ycoord.bias = c(-0.4,0.4),ycoord.var = c(8,13))

## 2D est2-dnorm ====
load("TMLE-est2-dnorm/result.Rdata")
p.2D.est2.dnorm <- plot.tmle(r'($\psi_2(\hat{Q}^*)$ - dnorm)',ycoord.var = c(8,13))

## 2D est3 ====
load("TMLE-est2b/result.Rdata")
p.2D.est3 <- plot.tmle(r'($\psi_{2b}(\hat{Q}^*)$)',ycoord.var = c(8.5,12.5))

p.2D <- plot_grid(
  p.2D.est2
  ,p.2D.est2.dnorm
  ,p.2D.est3
  , align = "hv"
  , ncol = 1
)

## 2D onestep-densratio ====
load("Onestep-est2a/result.Rdata")
p.2D.1densratio <- plot.tmle(r'($\psi_{2a}^{+}(\hat{Q})$)',ycoord.bias = c(-0.4,0.4),ycoord.var = c(8,13))

## 2D onestep-dnorm ====
load("Onestep-est2-dnorm/result.Rdata")
p.2D.1dnorm.sr <- plot.tmle(r'($\psi_2^{+}(\hat{Q})$ - dnorm)',ycoord.var = c(8,13))

## 2D onestep-bayes ====
load("Onestep-est2b/result.Rdata")
p.2D.1bayes <- plot.tmle(r'($\psi_{2b}^{+}(\hat{Q})$)',ycoord.var = c(8.5,12.5))

p.con1 <- plot_grid(
  p.2D.est2
  ,p.2D.est2.dnorm
  ,p.2D.est3
  , align = "hv"
  , ncol = 1
)

p.con2 <- plot_grid(
  p.2D.1densratio
  ,p.2D.1dnorm.sr
  ,p.2D.1bayes
  , align = "hv"
  , ncol = 1
)


p <- plot_grid(
  p.con1
  ,p.con2,
  align="hv",ncol=2)


ggsave("plot.pdf", plot = p, width = 16, height = 20, units = "in")



################################################ ATT ###################################################################
### multiM-d4 case ====
# the truth
load("../DGPs/1-truth-multi-d4-ATT.Rdata")

## 2D est2 ====
load("TMLE-est2a/ATT_result.Rdata")
p.2D.est2 <- plot.tmle(r'($\beta_{a}(\hat{Q}^*)$)')

## 2D est2-dnorm ====
load("TMLE-est2-dnorm/ATT_result.Rdata")
p.2D.est2.dnorm <- plot.tmle(r'($\beta(\hat{Q}^*)$ - dnorm)')

## 2D est3 ====
load("TMLE-est2b/ATT_result.Rdata")
p.2D.est3 <- plot.tmle(r'($\beta_{b}(\hat{Q}^*)$)')

p.2D <- plot_grid(
  p.2D.est2
  ,p.2D.est2.dnorm
  ,p.2D.est3
  , align = "hv"
  , ncol = 1
)

## 2D onestep-densratio ====
load("Onestep-est2a/ATT_result.Rdata")
p.2D.1densratio <- plot.tmle(r'($\beta_{a}^{+}(\hat{Q})$)')

## 2D onestep-dnorm ====
load("Onestep-est2-dnorm/ATT_result.Rdata")
p.2D.1dnorm.sr <- plot.tmle(r'($\beta^{+}(\hat{Q})$ - dnorm)')

## 2D onestep-bayes ====
load("Onestep-est2b/ATT_result.Rdata")
p.2D.1bayes <- plot.tmle(r'($\beta_{b}^{+}(\hat{Q})$)')

p.con1 <- plot_grid(
  p.2D.est2
  ,p.2D.est2.dnorm
  ,p.2D.est3
  , align = "hv"
  , ncol = 1
)

p.con2 <- plot_grid(
  p.2D.1densratio
  ,p.2D.1dnorm.sr
  ,p.2D.1bayes
  , align = "hv"
  , ncol = 1
)


p <- plot_grid(
  p.con1
  ,p.con2,
  align="hv",ncol=2)


ggsave("ATT_plot_d4.pdf", plot = p, width = 16, height = 20, units = "in")

