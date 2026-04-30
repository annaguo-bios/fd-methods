# sample size
n.vec <- c(500,1000,2000,4000,8000)

# number of simulations
nsim <- 1000



dgp.f.name="6-dgp-binaryM-continuousZ-interAZ.R" # name for the DGP function
truth = "6-truth-binaryM-continuousZ-tn06_07-interAZ.Rdata" # name for the truth.Rdata
out.path.TMLE= "TMLE/output/" # path for the output folder
out.path.Onestep= "Onestep/output/" # path for the output folder
pz <- "tn06_07"

# superlearner.L = "T"
# superlearner.M = "T"
# lib.L = "\"c('SL.glm','SL.bayesglm', 'SL.gam','SL.earth','SL.ranger','SL.svm','SL.xgboost','SL.mean')\"" # superlearner lib for L via bayes
# lib.M = "\"c('SL.glm','SL.bayesglm', 'SL.gam','SL.earth','SL.ranger','SL.svm','SL.xgboost','SL.mean')\"" # superlearner lib for M via bayes

for (i in seq_along(n.vec)){
  joblist <- c()
  for (t in 1:nsim){
    
    #if (i==1 & t %in% c(339)){t <- 2*nsim+t} # there are few outliers for n=250, change seed to avoid
    
    job <- paste0("Rscript main.R ",n.vec[i]," ",t," ", dgp.f.name," ", truth," ",out.path.TMLE," ", out.path.Onestep," ",pz)

    joblist <- c(joblist,job)
  }
  write.table(joblist, file = paste0("joblist_n",i,".txt") ,quote = F, col.names = F, row.names = F)
}
