args = commandArgs(trailingOnly=T)
m=as.integer(args[1]) # m-th multiple imputation

library(fdcausal)
library(SuperLearner)
library(dplyr)
library(mice)
library(boot)
# 
# dt1 <- read.csv("../../data/B_PROUD/export_mi.csv") %>% filter(.imp==0)
# dt2 <- read.csv("../../data/B_PROUD/export_continuousM.csv") %>% mutate(mds_alarm_to_needle = if_else(mds_tpa == 0, 0, mds_alarm_to_needle))
# dt2.imputed <- mice(dt2, m = 10, method = 'pmm', seed = 7)
# dt2.imputed_data <- complete(dt2.imputed, 'long')
# 
# dt <- dt2.imputed_data %>% 
#   # select(-X) %>% 
#   rename(A = STEMO, Y = mrs, M = mds_alarm_to_needle, X1 = systolic_before, X2 = nihss_before) %>% 
#   mutate(A = as.numeric(A), Y = as.numeric(Y)) %>% 
#   mutate(M1 = mds_tpa, M2=as.numeric(mds_tpa*M))
# 
# save(list = c('dt'), file="B_PROUD/dt.Rdata")


# bootstrap function for calculate DIM and LOR and the PMF of Y(a)
TMLE_boot <- function(dt,i){
  
  # resampling
  dt <- dt[i,]
  
  dt_Y0 <- dt %>% mutate(Y=as.numeric(Y==0))
  dt_Y1 <- dt %>% filter(Y>0) %>% mutate(Y=as.numeric(Y==1))
  dt_Y2 <- dt %>% filter(Y>1) %>% mutate(Y=as.numeric(Y==2))
  dt_Y3 <- dt %>% filter(Y>2) %>% mutate(Y=as.numeric(Y==3))
  dt_Y4 <- dt %>% filter(Y>3) %>% mutate(Y=as.numeric(Y==4))
  dt_Y5 <- dt %>% filter(Y>4) %>% mutate(Y=as.numeric(Y==5))
  
  ## get E[Y(a)=0]
  #a=1
  est_Y0_a1 <- estfd(a=1,data=dt_Y0,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  #a=0
  est_Y0_a0 <- estfd(a=0,data=dt_Y0,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y0_a1 <- est_Y0_a1$TMLE$estimated_psi
  tmle_Y0_a0 <- est_Y0_a0$TMLE$estimated_psi
  onestep_Y0_a1 <- est_Y0_a1$Onestep$estimated_psi
  onestep_Y0_a0 <- est_Y0_a0$Onestep$estimated_psi
  
  ## get E[Y(a)==1 | Y(a)>0]
  #a=1
  est_Y1_a1 <- estfd(a=1,data=dt_Y1,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y1_a0 <- estfd(a=0,data=dt_Y1,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y1_a1_conditional <- est_Y1_a1$TMLE$estimated_psi
  tmle_Y1_a0_conditional <- est_Y1_a0$TMLE$estimated_psi
  onestep_Y1_a1_conditional <- est_Y1_a1$Onestep$estimated_psi
  onestep_Y1_a0_conditional <- est_Y1_a0$Onestep$estimated_psi
  
  tmle_Y1_a1 <- tmle_Y1_a1_conditional * (1-tmle_Y0_a1)
  tmle_Y1_a0 <- tmle_Y1_a0_conditional * (1-tmle_Y0_a0)
  onestep_Y1_a1 <- onestep_Y1_a1_conditional * (1-onestep_Y0_a1)
  onestep_Y1_a0 <- onestep_Y1_a0_conditional * (1-onestep_Y0_a0)
  
  ## get E[Y(a)==2 | Y(a)>1]
  #a=1
  est_Y2_a1 <- estfd(a=1,data=dt_Y2,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y2_a0 <- estfd(a=0,data=dt_Y2,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y2_a1_conditional <- est_Y2_a1$TMLE$estimated_psi
  tmle_Y2_a0_conditional <- est_Y2_a0$TMLE$estimated_psi
  onestep_Y2_a1_conditional <- est_Y2_a1$Onestep$estimated_psi
  onestep_Y2_a0_conditional <- est_Y2_a0$Onestep$estimated_psi
  
  tmle_Y2_a1 <- tmle_Y2_a1_conditional * (1-tmle_Y1_a1 - tmle_Y0_a1)
  tmle_Y2_a0 <- tmle_Y2_a0_conditional * (1-tmle_Y1_a0 - tmle_Y0_a0)
  onestep_Y2_a1 <- onestep_Y2_a1_conditional * (1-onestep_Y1_a1 - onestep_Y0_a1)
  onestep_Y2_a0 <- onestep_Y2_a0_conditional * (1-onestep_Y1_a0 - onestep_Y0_a0)
  
  ## get E[Y(a)==3 | Y(a)>2]
  #a=1
  est_Y3_a1 <- estfd(a=1,data=dt_Y3,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y3_a0 <- estfd(a=0,data=dt_Y3,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y3_a1_conditional <- est_Y3_a1$TMLE$estimated_psi
  tmle_Y3_a0_conditional <- est_Y3_a0$TMLE$estimated_psi
  onestep_Y3_a1_conditional <- est_Y3_a1$Onestep$estimated_psi
  onestep_Y3_a0_conditional <- est_Y3_a0$Onestep$estimated_psi
  
  tmle_Y3_a1 <- tmle_Y3_a1_conditional * (1-tmle_Y2_a1 - tmle_Y1_a1 - tmle_Y0_a1)
  tmle_Y3_a0 <- tmle_Y3_a0_conditional * (1-tmle_Y2_a0 - tmle_Y1_a0 - tmle_Y0_a0)
  onestep_Y3_a1 <- onestep_Y3_a1_conditional * (1-onestep_Y2_a1 - onestep_Y1_a1 - onestep_Y0_a1)
  onestep_Y3_a0 <- onestep_Y3_a0_conditional * (1-onestep_Y2_a0 - onestep_Y1_a0 - onestep_Y0_a0)
  
  ## get E[Y(a)==4 | Y(a)>3]
  #a=1
  est_Y4_a1 <- estfd(a=1,data=dt_Y4,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y4_a0 <- estfd(a=0,data=dt_Y4,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y4_a1_conditional <- est_Y4_a1$TMLE$estimated_psi
  tmle_Y4_a0_conditional <- est_Y4_a0$TMLE$estimated_psi
  onestep_Y4_a1_conditional <- est_Y4_a1$Onestep$estimated_psi
  onestep_Y4_a0_conditional <- est_Y4_a0$Onestep$estimated_psi
  
  tmle_Y4_a1 <- tmle_Y4_a1_conditional * (1-tmle_Y3_a1 - tmle_Y2_a1 - tmle_Y1_a1 - tmle_Y0_a1)
  tmle_Y4_a0 <- tmle_Y4_a0_conditional * (1-tmle_Y3_a0 - tmle_Y2_a0 - tmle_Y1_a0 - tmle_Y0_a0)
  onestep_Y4_a1 <- onestep_Y4_a1_conditional * (1-onestep_Y3_a1 - onestep_Y2_a1 - onestep_Y1_a1 - onestep_Y0_a1)
  onestep_Y4_a0 <- onestep_Y4_a0_conditional * (1-onestep_Y3_a0 - onestep_Y2_a0 - onestep_Y1_a0 - onestep_Y0_a0)
  
  ## get E[Y(a)==5 | Y(a)>4]
  #a=1
  est_Y5_a1 <- estfd(a=1,data=dt_Y5,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y5_a0 <- estfd(a=0,data=dt_Y5,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), onestep=T, mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y5_a1_conditional <- est_Y5_a1$TMLE$estimated_psi
  tmle_Y5_a0_conditional <- est_Y5_a0$TMLE$estimated_psi
  onestep_Y5_a1_conditional <- est_Y5_a1$Onestep$estimated_psi
  onestep_Y5_a0_conditional <- est_Y5_a0$Onestep$estimated_psi
  
  tmle_Y5_a1 <- tmle_Y5_a1_conditional * (1-tmle_Y4_a1 - tmle_Y3_a1 - tmle_Y2_a1 - tmle_Y1_a1 - tmle_Y0_a1)
  tmle_Y5_a0 <- tmle_Y5_a0_conditional * (1-tmle_Y4_a0 - tmle_Y3_a0 - tmle_Y2_a0 - tmle_Y1_a0 - tmle_Y0_a0)
  onestep_Y5_a1 <- onestep_Y5_a1_conditional * (1-onestep_Y4_a1 - onestep_Y3_a1 - onestep_Y2_a1 - onestep_Y1_a1 - onestep_Y0_a1)
  onestep_Y5_a0 <- onestep_Y5_a0_conditional * (1-onestep_Y4_a0 - onestep_Y3_a0 - onestep_Y2_a0 - onestep_Y1_a0 - onestep_Y0_a0)
  
  ## get E[Y(a)==6]
  tmle_Y6_a1 <- 1 - tmle_Y5_a1 -tmle_Y4_a1 - tmle_Y3_a1 - tmle_Y2_a1 - tmle_Y1_a1 - tmle_Y0_a1
  tmle_Y6_a0 <- 1 - tmle_Y5_a0 -tmle_Y4_a0 - tmle_Y3_a0 - tmle_Y2_a0 - tmle_Y1_a0 - tmle_Y0_a0
  onestep_Y6_a1 <- 1 - onestep_Y5_a1 -onestep_Y4_a1 - onestep_Y3_a1 - onestep_Y2_a1 - onestep_Y1_a1 - onestep_Y0_a1
  onestep_Y6_a0 <- 1 - onestep_Y5_a0 -onestep_Y4_a0 - onestep_Y3_a0 - onestep_Y2_a0 - onestep_Y1_a0 - onestep_Y0_a0
  
  ## summarise the PMF
  tmle_Y_geq1_a1 <- tmle_Y1_a1 + tmle_Y2_a1 + tmle_Y3_a1 + tmle_Y4_a1 + tmle_Y5_a1 + tmle_Y6_a1 # p(Y(a=1)>=1)
  tmle_Y_geq1_a0 <- tmle_Y1_a0 + tmle_Y2_a0 + tmle_Y3_a0 + tmle_Y4_a0 + tmle_Y5_a0 + tmle_Y6_a0 # p(Y(a=0)>=1)
  onestep_Y_geq1_a1 <- onestep_Y1_a1 + onestep_Y2_a1 + onestep_Y3_a1 + onestep_Y4_a1 + onestep_Y5_a1 + onestep_Y6_a1 # p(Y(a=1)>=1)
  onestep_Y_geq1_a0 <- onestep_Y1_a0 + onestep_Y2_a0 + onestep_Y3_a0 + onestep_Y4_a0 + onestep_Y5_a0 + onestep_Y6_a0 # p(Y(a=0)>=1)
  
  tmle_Y_geq2_a1 <- tmle_Y2_a1 + tmle_Y3_a1 + tmle_Y4_a1 + tmle_Y5_a1 + tmle_Y6_a1 # p(Y(a=1)>=2)
  tmle_Y_geq2_a0 <- tmle_Y2_a0 + tmle_Y3_a0 + tmle_Y4_a0 + tmle_Y5_a0 + tmle_Y6_a0 # p(Y(a=0)>=2)
  onestep_Y_geq2_a1 <- onestep_Y2_a1 + onestep_Y3_a1 + onestep_Y4_a1 + onestep_Y5_a1 + onestep_Y6_a1 # p(Y(a=1)>=2)
  onestep_Y_geq2_a0 <- onestep_Y2_a0 + onestep_Y3_a0 + onestep_Y4_a0 + onestep_Y5_a0 + onestep_Y6_a0 # p(Y(a=0)>=2)
  
  tmle_Y_geq3_a1 <- tmle_Y3_a1 + tmle_Y4_a1 + tmle_Y5_a1 + tmle_Y6_a1 # p(Y(a=1)>=3)
  tmle_Y_geq3_a0 <- tmle_Y3_a0 + tmle_Y4_a0 + tmle_Y5_a0 + tmle_Y6_a0 # p(Y(a=0)>=3)
  onestep_Y_geq3_a1 <- onestep_Y3_a1 + onestep_Y4_a1 + onestep_Y5_a1 + onestep_Y6_a1 # p(Y(a=1)>=3)
  onestep_Y_geq3_a0 <- onestep_Y3_a0 + onestep_Y4_a0 + onestep_Y5_a0 + onestep_Y6_a0 # p(Y(a=0)>=3)
  
  tmle_Y_geq4_a1 <- tmle_Y4_a1 + tmle_Y5_a1 + tmle_Y6_a1 # p(Y(a=1)>=4)
  tmle_Y_geq4_a0 <- tmle_Y4_a0 + tmle_Y5_a0 + tmle_Y6_a0 # p(Y(a=0)>=4)
  onestep_Y_geq4_a1 <- onestep_Y4_a1 + onestep_Y5_a1 + onestep_Y6_a1 # p(Y(a=1)>=4)
  onestep_Y_geq4_a0 <- onestep_Y4_a0 + onestep_Y5_a0 + onestep_Y6_a0 # p(Y(a=0)>=4)
  
  tmle_Y_geq5_a1 <- tmle_Y5_a1 + tmle_Y6_a1 # p(Y(a=1)>=5)
  tmle_Y_geq5_a0 <- tmle_Y5_a0 + tmle_Y6_a0 # p(Y(a=0)>=5)
  onestep_Y_geq5_a1 <- onestep_Y5_a1 + onestep_Y6_a1 # p(Y(a=1)>=5)
  onestep_Y_geq5_a0 <- onestep_Y5_a0 + onestep_Y6_a0 # p(Y(a=0)>=5)
  
  tmle_Y_geq6_a1 <- tmle_Y6_a1 # p(Y(a=1)>=6)
  tmle_Y_geq6_a0 <- tmle_Y6_a0 # p(Y(a=0)>=6)
  onestep_Y_geq6_a1 <- onestep_Y6_a1 # p(Y(a=1)>=6)
  onestep_Y_geq6_a0 <- onestep_Y6_a0 # p(Y(a=0)>=6)
  
  ## target estimand
  # Difference in means (DIM)
  dim_tmle_a1 <- 0*tmle_Y0_a1 + 1*tmle_Y1_a1 + 2*tmle_Y2_a1 + 3*tmle_Y3_a1 + 4*tmle_Y4_a1 + 5*tmle_Y5_a1 + 6*tmle_Y6_a1
  dim_tmle_a0 <- 0*tmle_Y0_a0 + 1*tmle_Y1_a0 + 2*tmle_Y2_a0 + 3*tmle_Y3_a0 + 4*tmle_Y4_a0 + 5*tmle_Y5_a0 + 6*tmle_Y6_a0
  dim_tmle <- dim_tmle_a1 - dim_tmle_a0
  
  dim_onestep_a1 <- 0*onestep_Y0_a1 + 1*onestep_Y1_a1 + 2*onestep_Y2_a1 + 3*onestep_Y3_a1 + 4*onestep_Y4_a1 + 5*onestep_Y5_a1 + 6*onestep_Y6_a1
  dim_onestep_a0 <- 0*onestep_Y0_a0 + 1*onestep_Y1_a0 + 2*onestep_Y2_a0 + 3*onestep_Y3_a0 + 4*onestep_Y4_a0 + 5*onestep_Y5_a0 + 6*onestep_Y6_a0
  dim_onestep <- dim_onestep_a1 - dim_onestep_a0
  
  # Log-odds ratio (LOR)
  lor_tmle <- mean(c(
    log({tmle_Y_geq1_a1/(1-tmle_Y_geq1_a1)}/{tmle_Y_geq1_a0/(1-tmle_Y_geq1_a0)}),
    log({tmle_Y_geq2_a1/(1-tmle_Y_geq2_a1)}/{tmle_Y_geq2_a0/(1-tmle_Y_geq2_a0)}),
    log({tmle_Y_geq3_a1/(1-tmle_Y_geq3_a1)}/{tmle_Y_geq3_a0/(1-tmle_Y_geq3_a0)}),
    log({tmle_Y_geq4_a1/(1-tmle_Y_geq4_a1)}/{tmle_Y_geq4_a0/(1-tmle_Y_geq4_a0)}),
    log({tmle_Y_geq5_a1/(1-tmle_Y_geq5_a1)}/{tmle_Y_geq5_a0/(1-tmle_Y_geq5_a0)}),
    log({tmle_Y_geq6_a1/(1-tmle_Y_geq6_a1)}/{tmle_Y_geq6_a0/(1-tmle_Y_geq6_a0)})
  )
  )
  
  lor_onestep <- mean(c(
    log({onestep_Y_geq1_a1/(1-onestep_Y_geq1_a1)}/{onestep_Y_geq1_a0/(1-onestep_Y_geq1_a0)}),
    log({onestep_Y_geq2_a1/(1-onestep_Y_geq2_a1)}/{onestep_Y_geq2_a0/(1-onestep_Y_geq2_a0)}),
    log({onestep_Y_geq3_a1/(1-onestep_Y_geq3_a1)}/{onestep_Y_geq3_a0/(1-onestep_Y_geq3_a0)}),
    log({onestep_Y_geq4_a1/(1-onestep_Y_geq4_a1)}/{onestep_Y_geq4_a0/(1-onestep_Y_geq4_a0)}),
    log({onestep_Y_geq5_a1/(1-onestep_Y_geq5_a1)}/{onestep_Y_geq5_a0/(1-onestep_Y_geq5_a0)}),
    log({onestep_Y_geq6_a1/(1-onestep_Y_geq6_a1)}/{onestep_Y_geq6_a0/(1-onestep_Y_geq6_a0)})
  )
  )
  
  return(c(dim_tmle, dim_onestep, lor_tmle, lor_onestep,
           tmle_Y0_a1, tmle_Y1_a1, tmle_Y2_a1, tmle_Y3_a1, tmle_Y4_a1, tmle_Y5_a1, tmle_Y6_a1,
           tmle_Y0_a0, tmle_Y1_a0, tmle_Y2_a0, tmle_Y3_a0, tmle_Y4_a0, tmle_Y5_a0, tmle_Y6_a0,
           onestep_Y0_a1, onestep_Y1_a1, onestep_Y2_a1, onestep_Y3_a1, onestep_Y4_a1, onestep_Y5_a1, onestep_Y6_a1,
           onestep_Y0_a0, onestep_Y1_a0, onestep_Y2_a0, onestep_Y3_a0, onestep_Y4_a0, onestep_Y5_a0, onestep_Y6_a0))
  
  
  
}


set.seed(7)

load('dt.Rdata')

print('m: ');print(m)
dt.analysis <- dt %>% filter(.imp==m)
n <- nrow(dt.analysis)

# bootstrap for DIM and LOR
estimation.boot <- boot(dt.analysis, statistic=TMLE_boot, R=200)

set1.TMLE.dim.est <- mean(estimation.boot$t[,1])
set1.TMLE.lor.est <- mean(estimation.boot$t[,3])
set1.onestep.dim.est <- mean(estimation.boot$t[,2])
set1.onestep.lor.est <- mean(estimation.boot$t[,4])

set1.TMLE.dim.var <- var(estimation.boot$t[,1])
set1.TMLE.lor.var <- var(estimation.boot$t[,3])
set1.onestep.dim.var <- var(estimation.boot$t[,2])
set1.onestep.lor.var <- var(estimation.boot$t[,4])

save(list = c("estimation.boot", "set1.TMLE.dim.est", "set1.TMLE.lor.est", "set1.onestep.dim.est", "set1.onestep.lor.est",
              "set1.TMLE.dim.var", "set1.TMLE.lor.var", "set1.onestep.dim.var", "set1.onestep.lor.var"),
     file = paste0("/output/output_m",m,".Rdata"))
