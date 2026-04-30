args = commandArgs(trailingOnly=T)
out.path =args[1] # path for the output folder
out.name=args[2] # require specifying output file name, e.g. out.name <- "continuous_dat.Rdata"
n.vec.ind=args[3] # sample size indicator, 1: the long sample size vector, 0: the short sample size vector
nsim=as.integer(args[4]) # require specifying number of simulations, e.g. nsim <- 1000
truth=args[5] # path+name for the truth.Rdata
method=args[6] # method for estimation, e.g. method <- "Onestep/" or method <- "TMLE/"

# out.path ="output/"
# out.name= "result.Rdata"
# n.vec.ind=4
# nsim=1000
# truth="../DGPs/6-truth-binaryM-binaryZ.Rdata"
# method="TMLE/"

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)
library(stats)
library(xtable)
library(here)

if (n.vec.ind=="1"){
  n.vec <- c(250,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
} else if (n.vec.ind=="0"){
  n.vec <- c(250,500,1000,2000,4000)
} else if (n.vec.ind=="2"){
  n.vec <- c(500,1000,2000,4000)
} else if (n.vec.ind=="3"){
  n.vec <- c(500,1000,2000)
} else if (n.vec.ind=="4"){
  n.vec <- c(250,500,1000,2000,4000,8000)
} else if (n.vec.ind=="5"){
  n.vec <- c(500,1000,2000,4000,8000)
}else {
  n.vec <- c(as.integer(n.vec.ind))
}

#################################################
# Load truth
#################################################
load(truth)


############################
# Organize results
#
############################


# Load a random output file to get the levels of z
load(paste0(method, out.path, "output_",n.vec[1],"_1.Rdata"))


# record results
for (est in levels.z){
  assign(paste0('bias_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # bias
  assign(paste0('var_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # mse
  assign(paste0('est_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # estimate
  assign(paste0('est_matrix_EY1.', est), matrix(nrow = nsim, ncol = length(n.vec))) # estimate
  assign(paste0('est_matrix_EY0.', est), matrix(nrow = nsim, ncol = length(n.vec))) # estimate
  assign(paste0('ci_matrix_ate_lower.', est), matrix(nrow = nsim, ncol = length(n.vec))) # lower bound of 95% CI
  assign(paste0('ci_matrix_ate_upper.', est), matrix(nrow = nsim, ncol = length(n.vec))) # upper bound of 95% CI

  # 95% CI coverage
  assign(paste0('ci_coverage_ATE.', est), data.frame(n=n.vec, # sample size
                                                     coverage=vector(mode="integer",length=length(n.vec)), # 95% CI coverage
                                                     coverage.lower=vector(mode="integer",length=length(n.vec)),coverage.upper=vector(mode="integer",length=length(n.vec)))) # 95% CI for the above coverage
  # average point estimate
  assign(paste0('avg.ate.', est), data.frame(n=n.vec, # sample size
                                             est=vector(mode="integer",length=length(n.vec)), # E(phi_hat)
                                             sd=vector(mode="integer",length=length(n.vec)), # sd
                                             upper=vector(mode="integer",length=length(n.vec)), # E(phi_hat)-1.96*sd
                                             lower=vector(mode="integer",length=length(n.vec)))) # E(phi_hat)+1.96*sd)
  # average bias
  assign(paste0('avg.bias_ate.', est), data.frame(n=n.vec, # sample size
                                                  bias=vector(mode="integer",length=length(n.vec)), # sqrt(n)*E{E(phi_hat-phi_0)}
                                                  upper=vector(mode="integer",length=length(n.vec)), # bias-1.96*sd
                                                  lower=vector(mode="integer",length=length(n.vec)))) # bias+1.96*sd
  # average variance
  assign(paste0('avg.variance_ate.', est), data.frame(n=n.vec, # sample size
                                                      variance=vector(mode="integer",length=length(n.vec)), # n*E{E(phi_hat-E(phi_hat))^2}
                                                      upper=vector(mode="integer",length=length(n.vec)), # var-1.96*sd
                                                      lower=vector(mode="integer",length=length(n.vec)))) # var+1.96*sd

  # MSE
  assign(paste0('avg.MSE_ate.', est), data.frame(n=n.vec, # sample size
                                             mse=vector(mode="integer",length=length(n.vec)), # E{E(phi_hat-phi_0)^2}
                                             upper=vector(mode="integer",length=length(n.vec)), # mse-1.96*sd
                                             lower=vector(mode="integer",length=length(n.vec)))) # mse+1.96*sd

}



for (i in seq_along(n.vec)) {
  n <- n.vec[i]

  for (t in 1:nsim) {

    load(paste0(method, out.path, "output_", n, "_", t, ".Rdata"))

    for (est in levels.z) {

      # print(est)
      
      # Record bias
      bias_matrix <- get(paste0('bias_matrix_ate.', est))
      bias_matrix[t, i] <- get(paste0('hat_ATE.', est)) - ATE[est]
      assign(paste0('bias_matrix_ate.', est), bias_matrix, envir = parent.frame())
      
      # record variance
      var_matrix <- get(paste0('var_matrix_ate.', est))
      var_matrix[t, i] <- get(paste0('var_ATE.', est))
      assign(paste0('var_matrix_ate.', est), var_matrix, envir = parent.frame())

      # Record point estimate
      # ATE
      est_matrix <- get(paste0('est_matrix_ate.', est))
      est_matrix[t, i] <- get(paste0('hat_ATE.', est))
      assign(paste0('est_matrix_ate.', est), est_matrix, envir = parent.frame())
      
      if(est!='opt.linear'){
      # EY1
      est_matrix <- get(paste0('est_matrix_EY1.', est))
      est_matrix[t, i] <- get(paste0('hat_EY1.', est))
      assign(paste0('est_matrix_EY1.', est), est_matrix, envir = parent.frame())
      
      # EY0
      est_matrix <- get(paste0('est_matrix_EY0.', est))
      est_matrix[t, i] <- get(paste0('hat_EY0.', est))
      assign(paste0('est_matrix_EY0.', est), est_matrix, envir = parent.frame())
      }

      # Record lower CI
      ci_lower_matrix <- get(paste0('ci_matrix_ate_lower.', est))
      ci_lower_matrix[t, i] <- get(paste0('lower.ci_ATE.', est))
      assign(paste0('ci_matrix_ate_lower.', est), ci_lower_matrix, envir = parent.frame())

      # Record upper CI
      ci_upper_matrix <- get(paste0('ci_matrix_ate_upper.', est))
      ci_upper_matrix[t, i] <- get(paste0('upper.ci_ATE.', est))
      assign(paste0('ci_matrix_ate_upper.', est), ci_upper_matrix, envir = parent.frame())
    }
  }

  for (est in levels.z) {
    # Record CI coverage
    ci_coverage <- get(paste0('ci_coverage_ATE.', est))
    ci_coverage[i, "coverage"] <- mean((get(paste0('ci_matrix_ate_lower.', est))[, i] < ATE[est]) & (get(paste0('ci_matrix_ate_upper.', est))[, i] > ATE[est]))
    ci_coverage[i, "coverage.lower"] <- ci_coverage[i, "coverage"] - 1.96 * sqrt(ci_coverage[i, "coverage"] * (1 - ci_coverage[i, "coverage"]) / n)
    ci_coverage[i, "coverage.upper"] <- ci_coverage[i, "coverage"] + 1.96 * sqrt(ci_coverage[i, "coverage"] * (1 - ci_coverage[i, "coverage"]) / n)
    assign(paste0('ci_coverage_ATE.', est), ci_coverage, envir = parent.frame())

    # Record average point estimate
    avg_ate <- get(paste0('avg.ate.', est))
    avg_ate[i, "est"] <- mean(get(paste0('est_matrix_ate.', est))[, i])
    avg_ate[i, "sd"] <- sqrt(var(get(paste0('est_matrix_ate.', est))[, i]) / nsim)
    avg_ate[i, "upper"] <- avg_ate[i, "est"] + 1.96 * avg_ate[i, "sd"]
    avg_ate[i, "lower"] <- avg_ate[i, "est"] - 1.96 * avg_ate[i, "sd"]
    assign(paste0('avg.ate.', est), avg_ate, envir = parent.frame())

    # Record average bias
    avg_bias <- get(paste0('avg.bias_ate.', est))
    avg_bias[i, "bias"] <- sqrt(n) * mean(get(paste0('bias_matrix_ate.', est))[, i])
    avg_bias[i, "upper"] <- avg_bias[i, "bias"] + 1.96 * sqrt(var(sqrt(n) * get(paste0('bias_matrix_ate.', est))[, i]) / nsim)
    avg_bias[i, "lower"] <- avg_bias[i, "bias"] - 1.96 * sqrt(var(sqrt(n) * get(paste0('bias_matrix_ate.', est))[, i]) / nsim)
    assign(paste0('avg.bias_ate.', est), avg_bias, envir = parent.frame())

    # Record variance
    avg_variance <- get(paste0('avg.variance_ate.', est))
    avg_variance[i, "variance"] <- n * mean((get(paste0('bias_matrix_ate.', est))[, i] - mean(get(paste0('bias_matrix_ate.', est))[, i])) ^ 2)
    avg_variance[i, "upper"] <- avg_variance[i, "variance"] + 1.96 * sqrt(var(n * (get(paste0('bias_matrix_ate.', est))[, i] - mean(get(paste0('bias_matrix_ate.', est))[, i])) ^ 2) / nsim)
    avg_variance[i, "lower"] <- avg_variance[i, "variance"] - 1.96 * sqrt(var(n * (get(paste0('bias_matrix_ate.', est))[, i] - mean(get(paste0('bias_matrix_ate.', est))[, i])) ^ 2) / nsim)
    assign(paste0('avg.variance_ate.', est), avg_variance, envir = parent.frame())

    # Record MSE
    avg_mse <- get(paste0('avg.MSE_ate.', est))
    avg_mse[i, "mse"] <- n * mean(get(paste0('bias_matrix_ate.', est))[, i] ^ 2)
    avg_mse[i, "upper"] <- avg_mse[i, "mse"] + 1.96 * sqrt(var(n * (get(paste0('bias_matrix_ate.', est))[, i] ^ 2)) / nsim)
    avg_mse[i, "lower"] <- avg_mse[i, "mse"] - 1.96 * sqrt(var(n * (get(paste0('bias_matrix_ate.', est))[, i] ^ 2)) / nsim)
    assign(paste0('avg.MSE_ate.', est), avg_mse, envir = parent.frame())
  }

}

# Create output list
outlist <- unlist(lapply(levels.z, function(x) {
  paste0(c("bias_matrix_ate", "est_matrix_ate", "ci_matrix_ate_lower","var_matrix_ate","est_matrix_EY1", "est_matrix_EY0",
           "ci_matrix_ate_upper", "ci_coverage_ATE", "avg.ate",
           "avg.bias_ate", "avg.variance_ate", "avg.MSE_ate"),'.', x)
}))


# save data
save(list = outlist , file = paste0(method, out.name))

