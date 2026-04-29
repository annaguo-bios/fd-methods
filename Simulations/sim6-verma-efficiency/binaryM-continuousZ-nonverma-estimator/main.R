args = commandArgs(trailingOnly=T)
n=as.integer(args[1]) # sample size for the simulation
seed=as.integer(args[2])  # seeds for replication
dgp.f.name=args[3] # name for the DGP function
truth=args[4] # name for the truth.Rdata
out.path.tmle= args[5] # path for the output folder
out.path.onestep= args[6] # path for the output folder

# n=1000
# seed=1
# dgp.f.name="6-dgp-binaryM-binaryZ.R" # name for the DGP function
# truth = "6-truth-binaryM-binaryZ.Rdata" # name for the truth.Rdata
# out.path.tmle= "TMLE/output/" # path for the output folder
# out.path.onestep= "Onestep/output/" # path for the output folder


# by default the density ratio method is bayes

if(file.exists(paste0(out.path.tmle,"output_",n,"_",seed,".Rdata")) & file.exists(paste0(out.path.onestep,"output_",n,"_",seed,".Rdata")) ){
  cat("File exists!\n")
  quit(save = "no")  # Exit without error  
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
library(fdcausal)

set.seed(seed)

seed <- ifelse(seed%%1000==0, 1000, seed%%1000)
#################################################
# Functions
#################################################
source(paste0("../DGPs/",dgp.f.name)) # generate_data(n)

#################################################
# Load truth
#################################################
load(paste0("../DGPs/",truth))

# generate data
dat_output = generate_data(n)
data = dat_output$data
attach(data, warn.conflicts=FALSE)


# run TMLE
tmle_output_Y1 <- estfd(a=1,data=data,treatment='A', mediators='M', outcome='Y', covariates='Z',
                       estimator = c('onestep','tmle'), mediator.method='bayes', superlearner=F,crossfit=F,
                       linkA='identity')

print("Y1 done")
tmle_output_Y0 <- estfd(a=0,data=data,treatment='A', mediators='M', outcome='Y', covariates='Z',
                       estimator = c('onestep','tmle'), mediator.method='bayes', superlearner=F,crossfit=F,
                       linkA='identity')


# estimate E[Y(1)], E[Y(0)], and ATE
hat_E.Y1 = tmle_output_Y1$TMLE$estimated_psi
hat_E.Y0 = tmle_output_Y0$TMLE$estimated_psi
hat_ATE = hat_E.Y1 - hat_E.Y0

# lower CI
lower.ci_Y1 = tmle_output_Y1$TMLE$lower.ci
lower.ci_Y0 = tmle_output_Y0$TMLE$lower.ci
lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$TMLE$EIF-tmle_output_Y0$TMLE$EIF)^2)/n)

# upper CI
upper.ci_Y1 = tmle_output_Y1$TMLE$upper.ci
upper.ci_Y0 = tmle_output_Y0$TMLE$upper.ci
upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$TMLE$EIF-tmle_output_Y0$TMLE$EIF)^2)/n)

# compute bias
bias_Y1 = hat_E.Y1 - E.Y1[1]
bias_Y0 = hat_E.Y0 - E.Y0[1]
bias_ATE = hat_ATE - ATE[1]

save(list = c("tmle_output_Y1","tmle_output_Y0","bias_Y1","bias_Y0","bias_ATE","hat_E.Y1","hat_E.Y0","hat_ATE","lower.ci_Y1","lower.ci_Y0","lower.ci_ATE","upper.ci_Y1","upper.ci_Y0","upper.ci_ATE"),file = paste0(out.path.tmle,"output_",n,"_",seed,".Rdata"))



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

save(list = c("tmle_output_Y1","tmle_output_Y0","bias_Y1","bias_Y0","bias_ATE","hat_E.Y1","hat_E.Y0","hat_ATE","lower.ci_Y1","lower.ci_Y0","lower.ci_ATE","upper.ci_Y1","upper.ci_Y0","upper.ci_ATE"),file = paste0(out.path.onestep,"output_",n,"_",seed,".Rdata"))




