args = commandArgs(trailingOnly=T)
n=as.integer(args[1]) # sample size for the simulation
seed=as.integer(args[2])  # seeds for replication
dgp.f.name=args[3] # name for the DGP function
truth=args[4] # name for the truth.Rdata
out.path.TMLE= args[5] # path for the output folder
out.path.Onestep= args[6] # path for the output folder


# by default the density ratio method is bayes

if(file.exists(paste0(out.path.TMLE,"output_",n,"_",seed,".Rdata")) & file.exists(paste0(out.path.Onestep,"output_",n,"_",seed,".Rdata")) ){
  stop("File exists!")
}

library(ggplot2)
library(ggpubr)
library(latex2exp)
library(reshape2)
library(stats)
library(haldensify)
library(np)
library(xtable)
library(SuperLearner)
library(densratio)
library(MASS)
library(mvtnorm)
library(utils)

set.seed(seed)

seed <- ifelse(seed%%1000==0, 1000, seed%%1000)
#################################################
# Functions
#################################################
source(paste0("../DGPs/",dgp.f.name)) # generate_data(n)
source("../R/ATE.R")

#################################################
# Load truth
#################################################
load(paste0("../DGPs/",truth))

# generate data
dat_output = generate_data(n)
data = dat_output$data
attach(data, warn.conflicts=FALSE)


# run TMLE
output <- frontdoor_effect_est(c(1,0),data=data,treatment='A', mediator='M', outcome='Y', anchor='Z',
                                           onestep=TRUE, superlearner=TRUE,crossfit=FALSE,K=5,
                                           lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, cvg.criteria=0.01,
                                           formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                                           truncate_lower=0.01, truncate_upper=0.99,
                                           verbose=T)

print("Estimation done")

# levels of Z used for estimation
levels.z = unique(sub("^[^.]*\\.", "", names(output)[grep("^(TMLE|Onestep)\\.", names(output))]))

for (m in c('TMLE','Onestep')){ # loop over TMLE and one-step

  outlist <- c()

  for (est in levels.z){ # loop over levels of z

    out <- output[[paste0(m,'.',est)]] # get output for z

    assign(paste0('hat_ATE.',est), out$ATE) # estimated parameter
    assign(paste0('lower.ci_ATE.',est), out$lower.ci) # lower bound of 95% CI
    assign(paste0('upper.ci_ATE.',est), out$upper.ci) # upper bound of 95% CI
    assign(paste0('bias_ATE.',est), out$ATE - ATE[est]) # bias
    assign(paste0('var_ATE.',est), out$var_ATE) # bias

    outlist <- c(outlist, paste0('hat_ATE.',est), paste0('lower.ci_ATE.',est), paste0('upper.ci_ATE.',est), paste0('bias_ATE.',est), paste0('var_ATE.',est))
    
    out.Y1 <- output[["est.Y1"]][[m]][[paste0('out.',est)]] # get output for E(Y^1)
    out.Y0 <- output[["est.Y0"]][[m]][[paste0('out.',est)]] # get output for E(Y^0)
    
    assign(paste0('hat_EY1.',est), out.Y1$estimated) # estimated parameter
    assign(paste0('hat_EY0.',est), out.Y0$estimated) # estimated parameter
    
    outlist <- c(outlist, paste0('hat_EY1.',est), paste0('hat_EY0.',est))

  }

  # save the levels of Z such that it's easier in the organization of the output
  outlist <- c(outlist, "levels.z")

  save(list = outlist,file = paste0(get(paste0('out.path.',m)),"output_",n,"_",seed,".Rdata")) # save TMLE output

}



