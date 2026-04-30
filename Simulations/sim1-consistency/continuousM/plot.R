# packages
library(here)
setwd(here("./sim1-consistency/continuousM/"))

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
# the truth
load("../DGPs/1-truth-continuous.Rdata")

## continuous-est1 ====
load("TMLE-est1/result.Rdata")
p.con.est1 <- plot.tmle(r'($\psi_1(\hat{Q}^*)$)',ycoord.bias = c(-1.2,0.2))

## continuous-est1-dnorm ====
load("TMLE-est1-dnorm/result.Rdata")
p.con.est1.dnorm <- plot.tmle(r'($\psi_1(\hat{Q}^*)$ - dnorm)',ycoord.var = c(6.5,9.5))

## continuous-est2 ====
load("TMLE-est2a/result.Rdata")
p.con.est2 <- plot.tmle(r'($\psi_{2a}(\hat{Q}^*)$)',ycoord.bias = c(-0.2,0.3),ycoord.var = c(6,12.5))

## continuous-est2-dnorm ====
load("TMLE-est2-dnorm/result.Rdata")
p.con.est2.dnorm <- plot.tmle(r'($\psi_{2}(\hat{Q}^*)$ - dnorm)',ycoord.var = c(6.65,9.5))

## continuous-est3 ====
load("TMLE-est2b/result.Rdata")
p.con.est3 <- plot.tmle(r'($\psi_{2b}(\hat{Q}^*)$)',ycoord.var = c(6.4,9.6))

## continuous-onestep-np ====
load("Onestep-est1/result.Rdata")
p.con.1np <- plot.tmle(r'($\psi_1^{+}(\hat{Q})$)',ycoord.bias = c(-1.2,0.3))

## continuous-onestep-dnorm-sr ====
load("Onestep-est2-dnorm/result.Rdata")
p.con.1dnorm.sr <- plot.tmle(r'($\psi_{2}^{+}(\hat{Q})$ - dnorm)',ycoord.var=c(6.65,9.5))

## continuous-onestep-dnorm ====
load("Onestep-est1-dnorm/result.Rdata")
p.con.1dnorm <- plot.tmle(r'($\psi_1^{+}(\hat{Q})$ - dnorm)',ycoord.var = c(6.5,9.5))

## continuous-onestep-densratio ====
load("Onestep-est2a/result.Rdata")
p.con.1densratio <- plot.tmle(r'($\psi_{2a}^{+}(\hat{Q})$)',ycoord.bias = c(-0.2,0.3), ycoord.var = c(6,12.5))

## continuous-onestep-bayes ====
load("Onestep-est2b/result.Rdata")
p.con.1bayes <- plot.tmle(r'($\psi_{2b}^{+}(\hat{Q})$)',ycoord.var = c(6.4,9.6))

p.con1 <- plot_grid(
  p.con.est1
  ,p.con.est1.dnorm
  ,p.con.est2
  ,p.con.est2.dnorm
  ,p.con.est3
  , align = "hv"
  , ncol = 1
)

p.con2 <- plot_grid(
  p.con.1np
  ,p.con.1dnorm
  ,p.con.1densratio
  ,p.con.1dnorm.sr
  ,p.con.1bayes
  , align = "hv"
  , ncol = 1
)


p <- plot_grid(
  p.con1
  ,p.con2,
  align="hv",ncol=2)


ggsave("plot.pdf", plot = p, width = 16, height = 18, units = "in")

################################################ ATT ###################################################################
# the truth
load("../DGPs/1-truth-continuous-ATT.Rdata")

## continuous-est1 ====
load("TMLE-est1/ATT_result.Rdata")
p.con.est1 <- plot.tmle(r'($\beta_1(\hat{Q}^*)$)')

## continuous-est1-dnorm ====
load("TMLE-est1-dnorm/ATT_result.Rdata")
p.con.est1.dnorm <- plot.tmle(r'($\beta_1(\hat{Q}^*)$ - dnorm)')

## continuous-est2 ====
load("TMLE-est2a/ATT_result.Rdata")
p.con.est2 <- plot.tmle(r'($\beta_{a}(\hat{Q}^*)$)')

## continuous-est2-dnorm ====
load("TMLE-est2-dnorm/ATT_result.Rdata")
p.con.est2.dnorm <- plot.tmle(r'($\beta(\hat{Q}^*)$ - dnorm)')

## continuous-est3 ====
load("TMLE-est2b/ATT_result.Rdata")
p.con.est3 <- plot.tmle(r'($\beta_{b}(\hat{Q}^*)$)')

## continuous-onestep-np ====
load("Onestep-est1/ATT_result.Rdata")
p.con.1np <- plot.tmle(r'($\beta_1^{+}(\hat{Q})$)')

## continuous-onestep-dnorm-sr ====
load("Onestep-est2-dnorm/ATT_result.Rdata")
p.con.1dnorm.sr <- plot.tmle(r'($\psi^{+}(\hat{Q})$ - dnorm)')

## continuous-onestep-dnorm ====
load("Onestep-est1-dnorm/ATT_result.Rdata")
p.con.1dnorm <- plot.tmle(r'($\beta_1^{+}(\hat{Q})$ - dnorm)')

## continuous-onestep-densratio ====
load("Onestep-est2a/ATT_result.Rdata")
p.con.1densratio <- plot.tmle(r'($\beta_{a}^{+}(\hat{Q})$)')

## continuous-onestep-bayes ====
load("Onestep-est2b/ATT_result.Rdata")
p.con.1bayes <- plot.tmle(r'($\beta_{b}^{+}(\hat{Q})$)')

p.con1 <- plot_grid(
  p.con.est1
  ,p.con.est1.dnorm
  ,p.con.est2
  ,p.con.est2.dnorm
  ,p.con.est3
  , align = "hv"
  , ncol = 1
)

p.con2 <- plot_grid(
  p.con.1np
  ,p.con.1dnorm
  ,p.con.1densratio
  ,p.con.1dnorm.sr
  ,p.con.1bayes
  , align = "hv"
  , ncol = 1
)


p <- plot_grid(
  p.con1
  ,p.con2,
  align="hv",ncol=2)


ggsave("ATT_plot_continuous.pdf", plot = p, width = 16, height = 18, units = "in")
