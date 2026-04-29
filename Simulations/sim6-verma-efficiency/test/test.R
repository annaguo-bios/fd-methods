set.seed(7)

generate_data <- function(n,parA = c(0.3,0.2), parU=c(1,1,-1,0), parM = c(-1,1,1,2), parY = c(1, 1, 1, 0), sd.U=1, sd.Y=1){ # change the parM to c(-1,1,1,0)
  
  Z <- rnorm(n,1,1) # p(Z)
  
  A <- rbinom(n, 1, plogis(parA[1] + parA[2]*Z)) # p(A|X)
  
  U <- parU[1] + parU[2]*A + parU[3]*plogis(parA[1] + parA[2]*Z) + parU[4]*A*Z + rnorm(n,0,sd.U) # p(U|A,X)
  
  M <- rbinom(n,1,plogis(parM[1] + parM[2]*A + parM[3]*Z + parM[4]*A*Z)) # p(M|A,X)
  
  Y <- parY[1]*U + parY[2]*M + rnorm(n, 0, sd.Y) # p(Y|U,M)
  
  data <- data.frame(Z=Z, U=U, A=A, M=M, Y=Y)
  
  return(list(data = data, 
              parA=parA, 
              parU=parU,
              parM=parM, 
              parY=parY,
              sd.U=sd.U,
              sd.Y=sd.Y))
}

setwd(here::here('sim6-verma-efficiency/test'))
source("../R/ATE-continuousZ-binaryM.R")

n.vec = c(1000,2000,4000,8000) # sample size
nsim = 200 # number of simulations
pz='nonverma'

if (pz=='dgp-pz'){
  
  tilde_pz <- function(z) {dnorm(z,1,1)}
  
}else if (pz=='uniform100-pz'){
  
  tilde_pz <- function(z) { dunif(z,-100,100)}
  
}else if (pz=='normal01-pz'){
  
  tilde_pz <- function(z) { dnorm(z,0,1)}
  
}else if (pz=='normal10.1-pz'){
  
  tilde_pz <- function(z) {dnorm(z,10,1)}
  
}else if (pz=='mix0.9'){
  
  eps <- 0.9
  
  tilde_pz <- function(z) {eps*dnorm(z,10,1)+eps*dnorm(z,0,1)}
  
}else if (pz=='mix0.7'){
  
  eps <- 0.7
  
  tilde_pz <- function(z) {eps*dnorm(z,10,1)+eps*dnorm(z,0,1)}
  
}else if (pz=='tn06_07'){
  
  # keep sample_pz in the bulk range of N(1,1)
  a <- qnorm(0.001, 1, 1)
  b <- qnorm(0.999, 1, 1)
  
  mu_t <- 0.6
  sd_t <- 0.7
  
  Ct <- pnorm(b, mu_t, sd_t) - pnorm(a, mu_t, sd_t)
  
  tilde_pz <- function(z){
    dnorm(z, mu_t, sd_t) * (z >= a & z <= b) / Ct
  }
  
  # sample_pz by rejection sampling (guaranteed in [a,b])
  sample_pz <- numeric(0)
  while(length(sample_pz) < n_sample){
    z <- rnorm(n_sample, mu_t, sd_t)
    z <- z[z >= a & z <= b]
    sample_pz <- c(sample_pz, z)
  }
  sample_pz <- sample_pz[1:n_sample]
}

load(paste0('../test/6-truth-binaryM-continuousZ-',pz,'-local.Rdata')) # load true ATE for each z

for (est in levels.z){
  assign(paste0('bias_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # bias
  assign(paste0('var_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # mse
  assign(paste0('est_matrix_ate.', est), matrix(nrow = nsim, ncol = length(n.vec))) # estimate
  assign(paste0('ci_matrix_ate_lower.', est), matrix(nrow = nsim, ncol = length(n.vec))) # lower bound of 95% CI
  assign(paste0('ci_matrix_ate_upper.', est), matrix(nrow = nsim, ncol = length(n.vec))) # upper bound of 95% CI
  
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

for (i in seq_along(n.vec)){
  
  n = n.vec[i]
  
  for (t in 1:nsim){
    
    dat_output = generate_data(n)
    data = dat_output$data
    attach(data, warn.conflicts=FALSE)
    
    if(pz!='nonverma'){
      invisible(
        capture.output(output <- frontdoor_effect_est(c(1,0),data=data,treatment='A', mediator='M', outcome='Y', anchor='Z',tilde_pz=tilde_pz,
                                                      onestep=TRUE, superlearner=TRUE,crossfit=FALSE,K=5,
                                                      lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, cvg.criteria=0.01,
                                                      formulaY="Y ~ 1 + A + M + I(plogis(0.3 + 0.2 * Z))", formulaA="A ~ .", formulaM="M~A+Z+A*Z", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                                                      truncate_lower=0.01, truncate_upper=0.99,
                                                      verbose=T)
        ))
      
      if(t %% 50 == 0) {
        cat("Completed", t, "simulations for n =", n, "\n")
      }
      
      
      m = 'Onestep'
      est = levels.z
      
      out <- output[[paste0(m,'.',est)]] # get output for z
      
      assign(paste0('hat_ATE.',est), out$ATE) # estimated parameter
      assign(paste0('lower.ci_ATE.',est), out$lower.ci) # lower bound of 95% CI
      assign(paste0('upper.ci_ATE.',est), out$upper.ci) # upper bound of 95% CI
      assign(paste0('bias_ATE.',est), out$ATE - ATE) # bias
      assign(paste0('var_ATE.',est), out$var_ATE) # bias
      
      out.Y1 <- output[["est.Y1"]][[m]][[paste0('out.',est)]] # get output for E(Y^1)
      out.Y0 <- output[["est.Y0"]][[m]][[paste0('out.',est)]] # get output for E(Y^0)
      
      assign(paste0('hat_EY1.',est), out.Y1$estimated) # estimated parameter
      assign(paste0('hat_EY0.',est), out.Y0$estimated) # estimated parameter
      
      ## organize ##
      
      # Record bias
      bias_matrix <- get(paste0('bias_matrix_ate.', est))
      bias_matrix[t, i] <- get(paste0('hat_ATE.', est)) - ATE
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
      
      # Record lower CI
      ci_lower_matrix <- get(paste0('ci_matrix_ate_lower.', est))
      ci_lower_matrix[t, i] <- get(paste0('lower.ci_ATE.', est))
      assign(paste0('ci_matrix_ate_lower.', est), ci_lower_matrix, envir = parent.frame())
      
      # Record upper CI
      ci_upper_matrix <- get(paste0('ci_matrix_ate_upper.', est))
      ci_upper_matrix[t, i] <- get(paste0('upper.ci_ATE.', est))
      assign(paste0('ci_matrix_ate_upper.', est), ci_upper_matrix, envir = parent.frame())
    }else{
      # run TMLE
      invisible(
        capture.output(tmle_output_Y1 <- estfd(a=1,data=data,treatment='A', mediators='M', outcome='Y', covariates='Z',
                             estimator = c('onestep','tmle'), mediator.method='bayes', superlearner=F,crossfit=F,formulaM = "M~1+A+Z+A*Z")
        ))
      

      invisible(
        capture.output(tmle_output_Y0 <- estfd(a=0,data=data,treatment='A', mediators='M', outcome='Y', covariates='Z',
                             estimator = c('onestep','tmle'), mediator.method='bayes', superlearner=F,crossfit=F,formulaM = "M~1+A+Z+A*Z")
        ))
      
      if(t %% 50 == 0) {
        cat("Completed", t, "simulations for n =", n, "\n")
      }
      
      # estimate E[Y(1)], E[Y(0)], and ATE
      hat_E.Y1 = tmle_output_Y1$Onestep$estimated_psi
      hat_E.Y0 = tmle_output_Y0$Onestep$estimated_psi
      hat_ATE = hat_E.Y1 - hat_E.Y0
      
      # lower CI
      lower.ci_Y1 = tmle_output_Y1$Onestep$lower.ci
      lower.ci_Y0 = tmle_output_Y0$Onestep$lower.ci
      lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$Onestep$EIF-tmle_output_Y0$Onestep$EIF)^2)/n)
      
      # upper CI
      upper.ci_Y1 = tmle_output_Y1$Onestep$upper.ci
      upper.ci_Y0 = tmle_output_Y0$Onestep$upper.ci
      upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$Onestep$EIF-tmle_output_Y0$Onestep$EIF)^2)/n)
      
      # compute bias
      bias_Y1 = hat_E.Y1 - E.Y1[1]
      bias_Y0 = hat_E.Y0 - E.Y0[1]
      bias_ATE = hat_ATE - ATE[1]
      
      # compute variance
      var_ATE = mean((tmle_output_Y1$Onestep$EIF-tmle_output_Y0$Onestep$EIF)^2)
      
      
      ## organize ##
      
      # Record bias
      bias_matrix <- get(paste0('bias_matrix_ate.', est))
      bias_matrix[t, i] <- hat_ATE - ATE
      assign(paste0('bias_matrix_ate.', est), bias_matrix, envir = parent.frame())
      
      # record variance
      var_matrix <- get(paste0('var_matrix_ate.', est))
      var_matrix[t, i] <- var_ATE
      assign(paste0('var_matrix_ate.', est), var_matrix, envir = parent.frame())
      
      # Record point estimate
      # ATE
      est_matrix <- get(paste0('est_matrix_ate.', est))
      est_matrix[t, i] <- hat_ATE
      assign(paste0('est_matrix_ate.', est), est_matrix, envir = parent.frame())
      
      # Record lower CI
      ci_lower_matrix <- get(paste0('ci_matrix_ate_lower.', est))
      ci_lower_matrix[t, i] <- lower.ci_ATE
      assign(paste0('ci_matrix_ate_lower.', est), ci_lower_matrix, envir = parent.frame())
      
      # Record upper CI
      ci_upper_matrix <- get(paste0('ci_matrix_ate_upper.', est))
      ci_upper_matrix[t, i] <- upper.ci_ATE
      assign(paste0('ci_matrix_ate_upper.', est), ci_upper_matrix, envir = parent.frame())
      
      
    }
  }

  
}


for(i in seq_along(n.vec)) {
  n = n.vec[i]
  for (est in levels.z) {
    
    # # Record average point estimate
    # avg_ate <- get(paste0('avg.ate.', est))
    # avg_ate[i, "est"] <- mean(get(paste0('est_matrix_ate.', est))[, i])
    # avg_ate[i, "sd"] <- sqrt(var(get(paste0('est_matrix_ate.', est))[, i]) / nsim)
    # avg_ate[i, "upper"] <- avg_ate[i, "est"] + 1.96 * avg_ate[i, "sd"]
    # avg_ate[i, "lower"] <- avg_ate[i, "est"] - 1.96 * avg_ate[i, "sd"]
    # assign(paste0('avg.ate.', est), avg_ate, envir = parent.frame())
    
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

outlist <- unlist(lapply(levels.z, function(x) {
  paste0(c("bias_matrix_ate", "est_matrix_ate", "ci_matrix_ate_lower","var_matrix_ate",
           "ci_matrix_ate_upper", 
           "avg.bias_ate", "avg.variance_ate", "avg.MSE_ate"),'.', x)
}))


# save data
save(list = outlist , file = paste0('/Users/apple/Library/CloudStorage/Dropbox/Front-door_Anna/fd-methods/Simulations/sim6-verma-efficiency/test/',pz,'.Rdata'))
