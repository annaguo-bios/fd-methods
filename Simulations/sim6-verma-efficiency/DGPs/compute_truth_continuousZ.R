args = commandArgs(trailingOnly=T)
truth.f.name=args[1] # require specifying truth function, e.g. truth.f.name <- "6-truth-binaryM-continuousZ.R"
dgp.f.name=args[2] # the dgp.R function, e.g. dgp.f.name <- "6-dgp-binaryM-continuousZ-interAZ.R"
out.name=args[3] # require specifying the name for the output file, e.g. out.name <- "0-truth-binary.Rdata"
N=as.integer(args[4]) # require specifying the sample size for computing the truth numerically, e.g. N <- 100000
pz = args[5] # pz='dgp-pz'


tilde_pz <- NULL
n_sample <- 500

if (pz=='dgp-pz'){
  
  tilde_pz <- function(z) {dnorm(z,1,1)}
  
  sample_pz <- rnorm(n_sample,1,1)
  
}else if (pz=='uniform100-pz'){
  
  tilde_pz <- function(z) { dunif(z,-100,100)}
  
  sample_pz <- runif(n_sample,-100,100)
  
}else if (pz=='uniform1-pz'){
  
  tilde_pz <- function(z) { dunif(z,9,10)}
  
  sample_pz <- runif(n_sample,9,10)
  
}else if (pz=='normal01-pz'){
  
  tilde_pz <- function(z) { dnorm(z,0,1)}
  
  sample_pz <- rnorm(n_sample,0,1)
  
}else if (pz=='normal10.1-pz'){
  
  tilde_pz <- function(z) {dnorm(z,10,1)}
  
  sample_pz <- rnorm(n_sample,10,1)
  
}else if (pz=='mix'){
  
  eps <- 0.3  # try 0.01, 0.02, 0.05
  
  tilde_pz <- function(z) (1-eps)*dnorm(z,1,1) + eps*dnorm(z,10,1)
  
  mix_ind <- rbinom(n_sample, 1, eps)
  sample_pz <- ifelse(mix_ind==1, rnorm(n_sample,10,1), rnorm(n_sample,1,1))
  
}

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)
library(stats)
library(xtable)
library(cubature)
library(MASS)
library(mvtnorm)
library(Matrix)

set.seed(7)

#################################################
# Functions
#################################################
source(truth.f.name) # compute_truth(n)
source(dgp.f.name) # generate_data(n)


#################################################
# True parameters
#################################################

## true psi and variance
truth_output <- compute_truth(n = N, tilde_pz, sample_pz)

# to check that the mean of EIF=0
mean.EIF.Y1<- mean(truth_output$EIF_Y1); mean.EIF.Y1
mean.EIF.Y0<- mean(truth_output$EIF_Y0); mean.EIF.Y0

E.Y1 = truth_output$E.Y1; names(E.Y1) = c("all.z")
VAR.Y1 = truth_output$VAR.Y1; names(VAR.Y1) = c("all.z")
#
E.Y0 = truth_output$E.Y0; names(E.Y0) = c("all.z")
VAR.Y0 = truth_output$VAR.Y0; names(VAR.Y0) = c("all.z")
#
ATE = truth_output$ATE; names(ATE) = c("all.z")
VAR.ATE = truth_output$VAR.ATE; names(VAR.ATE) = c("all.z")

# print out the results
at.Y1 <- data.frame(E.Y1=E.Y1, VAR.Y1=VAR.Y1)
rownames(at.Y1) <- c("all.z")

at.Y0 <- data.frame(E.Y0=E.Y0, VAR.Y0=VAR.Y0)
rownames(at.Y0) <- c("all.z")

at.ATE <- data.frame(ATE=ATE, VAR.ATE=VAR.ATE)
rownames(at.ATE) <- c("all.z")


# print("E[Y1] and it's variance")
# print(at.Y1)
# 
# 
# print("E[Y0] and it's variance")
# print(at.Y0)

cat('tilde_pz: ', pz, '\n')
print("ATE and it's variance")
print(at.ATE)
cat('\n')


save(list = c("E.Y1","E.Y0","ATE","VAR.Y1","VAR.Y0","VAR.ATE","mean.EIF.Y1","mean.EIF.Y0"),file = out.name)
