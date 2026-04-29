args = commandArgs(trailingOnly=T)
n=as.integer(args[1]) # sample size for the simulation
seed=as.integer(args[2])  # seeds for replication
dgp.f.name=args[3] # name for the DGP function
truth=args[4] # name for the truth.Rdata
out.path.tmle= args[5] # path for the output folder
out.path.onestep= args[6] # path for the output folder
mediator.method=args[7]
superlearner=as.logical(args[8])
crossfit=as.logical(args[9])
K=as.integer(args[10])
treatment=args[11]
mediators=args[12]
outcome=args[13]
covariates=args[14]
lib=args[15]
linkA=args[16]

# Set optional argument with default value
ATT.arg = ifelse(length(args) >= 17, as.logical(args[17]), FALSE)

# Example: Print values
cat("ATT Argument (default FALSE):", ATT.arg, "\n")

# prefix depending on ATT
if (ATT.arg){prefix = "ATT_"} else {prefix = ""}


if(file.exists(paste0(out.path.tmle,prefix,"output_",n,"_",seed,".Rdata")) & file.exists(paste0(out.path.onestep,prefix,"output_",n,"_",seed,".Rdata"))){
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
library(fdcausal)


set.seed(seed)
#################################################
# Functions
#################################################
source(dgp.f.name) # generate_data(n)

#################################################
# Load truth
#################################################
load(truth)

# generate data
dat_output = generate_data(n)
data = dat_output$data
attach(data, warn.conflicts=FALSE)

# run TMLE
if (ATT.arg){
  
  Y <- data[[outcome]]
  A <- data[[treatment]]
  estimated_psi <- mean((A==1)/mean(A==1)*Y)
  EIF <- (A==1)/mean(A==1)*(Y- estimated_psi)
  
  lower.ci <- estimated_psi - 1.96*sqrt(mean(EIF^2)/n)
  upper.ci <- estimated_psi + 1.96*sqrt(mean(EIF^2)/n)
  
  Onestep <- list(estimated_psi=estimated_psi, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)
  TMLE <- list(estimated_psi=estimated_psi, EIF=EIF, lower.ci=lower.ci, upper.ci=upper.ci)
  
  tmle_output_Y1.F <- list(TMLE=TMLE, Onestep=Onestep)
  tmle_output_Y1.T <- list(TMLE=TMLE, Onestep=Onestep)
  
}else{
  
  # run TMLE
  tmle_output_Y1.F <- estfd(a=1,data=data,treatment=treatment, mediators=eval(parse(text = mediators)), outcome=outcome, covariates=eval(parse(text = covariates)),
                         mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                         lib = eval(parse(text = lib)), n.iter=15000, eps=F, cvg.criteria=n^{-1/2},
                         estimator='tmle',boundedsubmodelY=F)
  
  tmle_output_Y1.T <- estfd(a=1,data=data,treatment=treatment, mediators=eval(parse(text = mediators)), outcome=outcome, covariates=eval(parse(text = covariates)),
                           mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                           lib = eval(parse(text = lib)), n.iter=15000, eps=F, cvg.criteria=n^{-1/2},
                           estimator='tmle',boundedsubmodelY=T)
}



print("Y1 done")
tmle_output_Y0.F <- estfd(a=0,data=data,treatment=treatment, mediators=eval(parse(text = mediators)), outcome=outcome, covariates=eval(parse(text = covariates)),
                         mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                         lib = eval(parse(text = lib)), n.iter=15000, eps=F, cvg.criteria=n^{-1/2},
                         linkA=linkA,ATT=ATT.arg,
                         estimator='tmle',boundedsubmodelY=F)

tmle_output_Y0.T <- estfd(a=0,data=data,treatment=treatment, mediators=eval(parse(text = mediators)), outcome=outcome, covariates=eval(parse(text = covariates)),
                       mediator.method=mediator.method, superlearner=superlearner,crossfit=crossfit,K=K,
                       lib = eval(parse(text = lib)), n.iter=15000, eps=F, cvg.criteria=n^{-1/2},
                       linkA=linkA,ATT=ATT.arg,
                       estimator='tmle',boundedsubmodelY=T)

## without bounded submodel ==

# estimate E[Y(1)], E[Y(0)], and ATE
hat_E.Y1 = tmle_output_Y1.F$TMLE$estimated_psi
hat_E.Y0 = tmle_output_Y0.F$TMLE$estimated_psi
hat_ATE = hat_E.Y1 - hat_E.Y0

# lower CI
lower.ci_Y1 = tmle_output_Y1.F$TMLE$lower.ci
lower.ci_Y0 = tmle_output_Y0.F$TMLE$lower.ci
lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1.F$TMLE$EIF-tmle_output_Y0.F$TMLE$EIF)^2)/n)
# lower.ci_Y1 = hat_E.Y1-1.96*sqrt(mean(tmle_output_Y1$EIF^2-tmle_output_Y1$EIF)/n)
# lower.ci_Y0 = hat_E.Y0-1.96*sqrt(mean(tmle_output_Y0$EIF^2-tmle_output_Y0$EIF)/n)
# lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2-(tmle_output_Y1$EIF-tmle_output_Y0$EIF))/n)

# upper CI
upper.ci_Y1 = tmle_output_Y1.F$TMLE$upper.ci
upper.ci_Y0 = tmle_output_Y0.F$TMLE$upper.ci
upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1.F$TMLE$EIF-tmle_output_Y0.F$TMLE$EIF)^2)/n)
# lower.ci_Y1 = hat_E.Y1+1.96*sqrt(mean(tmle_output_Y1$EIF^2-tmle_output_Y1$EIF)/n)
# lower.ci_Y0 = hat_E.Y0+1.96*sqrt(mean(tmle_output_Y0$EIF^2-tmle_output_Y0$EIF)/n)
# upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2-(tmle_output_Y1$EIF-tmle_output_Y0$EIF))/n)

# compute bias
bias_Y1 = hat_E.Y1 - E.Y1
bias_Y0 = hat_E.Y0 - E.Y0
bias_ATE = hat_ATE - ATE

save(list = c("tmle_output_Y1.F","tmle_output_Y0.F","bias_Y1","bias_Y0","bias_ATE","hat_E.Y1","hat_E.Y0","hat_ATE","lower.ci_Y1","lower.ci_Y0","lower.ci_ATE","upper.ci_Y1","upper.ci_Y0","upper.ci_ATE"),file = paste0(out.path.tmle,prefix,"output_",n,"_",seed,".Rdata"))

## with bounded submodel ====
# estimate E[Y(1)], E[Y(0)], and ATE
hat_E.Y1 = tmle_output_Y1.T$TMLE$estimated_psi
hat_E.Y0 = tmle_output_Y0.T$TMLE$estimated_psi
hat_ATE = hat_E.Y1 - hat_E.Y0

# lower CI
lower.ci_Y1 = tmle_output_Y1.T$TMLE$lower.ci
lower.ci_Y0 = tmle_output_Y0.T$TMLE$lower.ci
lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1.T$TMLE$EIF-tmle_output_Y0.T$TMLE$EIF)^2)/n)
# lower.ci_Y1 = hat_E.Y1-1.96*sqrt(mean(tmle_output_Y1$EIF^2-tmle_output_Y1$EIF)/n)
# lower.ci_Y0 = hat_E.Y0-1.96*sqrt(mean(tmle_output_Y0$EIF^2-tmle_output_Y0$EIF)/n)
# lower.ci_ATE = hat_ATE - 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2-(tmle_output_Y1$EIF-tmle_output_Y0$EIF))/n)

# upper CI
upper.ci_Y1 = tmle_output_Y1.T$TMLE$upper.ci
upper.ci_Y0 = tmle_output_Y0.T$TMLE$upper.ci
upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1.T$TMLE$EIF-tmle_output_Y0.T$TMLE$EIF)^2)/n)
# lower.ci_Y1 = hat_E.Y1+1.96*sqrt(mean(tmle_output_Y1$EIF^2-tmle_output_Y1$EIF)/n)
# lower.ci_Y0 = hat_E.Y0+1.96*sqrt(mean(tmle_output_Y0$EIF^2-tmle_output_Y0$EIF)/n)
# upper.ci_ATE = hat_ATE + 1.96*sqrt(mean((tmle_output_Y1$EIF-tmle_output_Y0$EIF)^2-(tmle_output_Y1$EIF-tmle_output_Y0$EIF))/n)

# compute bias
bias_Y1 = hat_E.Y1 - E.Y1
bias_Y0 = hat_E.Y0 - E.Y0
bias_ATE = hat_ATE - ATE

save(list = c("tmle_output_Y1.T","tmle_output_Y0.T","bias_Y1","bias_Y0","bias_ATE","hat_E.Y1","hat_E.Y0","hat_ATE","lower.ci_Y1","lower.ci_Y0","lower.ci_ATE","upper.ci_Y1","upper.ci_Y0","upper.ci_ATE"),file = paste0(out.path.onestep,prefix,"output_",n,"_",seed,".Rdata"))
