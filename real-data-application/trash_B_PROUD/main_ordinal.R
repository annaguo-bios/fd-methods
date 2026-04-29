args = commandArgs(trailingOnly=T)
m=as.integer(args[1]) # m-th multiple imputation

library(fdcausal)
library(SuperLearner)
library(dplyr)
library(mice)
library(boot)

set.seed(7)

## data
dt2 <- read.csv("export_continuousM.csv") %>% mutate(mds_alarm_to_needle = ifelse(mds_tpa == 0, 0, mds_alarm_to_needle))
dt2.imputed <- mice(dt2, m = 10, method = 'pmm', seed = 7)
dt2.imputed_data <- complete(dt2.imputed, 'long')

dt <- dt2.imputed_data %>% 
  # select(-X) %>% 
  rename(A = STEMO, Y = mrs, M = mds_alarm_to_needle, X1 = systolic_before, X2 = nihss_before) %>% 
  mutate(A = as.numeric(A), Y = as.numeric(Y)) %>% 
  mutate(M1 = mds_tpa, M2=as.numeric(mds_tpa*M))



TMLE_m <- function(dt){
  
  n <- nrow(dt)
  
  dt_Y0 <- dt %>% mutate(Y=as.numeric(Y==0))
  dt_Y1 <- dt %>% mutate(Y=as.numeric(Y<=1))
  dt_Y2 <- dt %>% mutate(Y=as.numeric(Y<=2))
  dt_Y3 <- dt %>% mutate(Y=as.numeric(Y<=3))
  dt_Y4 <- dt %>% mutate(Y=as.numeric(Y<=4))
  dt_Y5 <- dt %>% mutate(Y=as.numeric(Y<=5))
  
  ## get E[Y(a)=0]
  #a=1
  est_Y0_a1 <- estfd(a=1,data=dt_Y0,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  #a=0
  est_Y0_a0 <- estfd(a=0,data=dt_Y0,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y0_a1 <- est_Y0_a1$TMLE$estimated_psi
  tmle_Y0_a0 <- est_Y0_a0$TMLE$estimated_psi
  tmle_Y0_a1_EIF <- est_Y0_a0$TMLE$EIF
  tmle_Y0_a0_EIF <- est_Y0_a0$TMLE$EIF
  onestep_Y0_a1 <- est_Y0_a1$Onestep$estimated_psi
  onestep_Y0_a0 <- est_Y0_a0$Onestep$estimated_psi
  onestep_Y0_a1_EIF <- est_Y0_a0$Onestep$EIF
  onestep_Y0_a0_EIF <- est_Y0_a0$Onestep$EIF
  
  ## get E[Y(a)<=1]
  #a=1
  est_Y1_a1 <- estfd(a=1,data=dt_Y1,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y1_a0 <- estfd(a=0,data=dt_Y1,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y1_a1_leq <- est_Y1_a1$TMLE$estimated_psi
  tmle_Y1_a0_leq <- est_Y1_a0$TMLE$estimated_psi
  tmle_Y1_a1_leq_EIF <- est_Y1_a1$TMLE$EIF
  tmle_Y1_a0_leq_EIF <- est_Y1_a0$TMLE$EIF
  onestep_Y1_a1_leq <- est_Y1_a1$Onestep$estimated_psi
  onestep_Y1_a0_leq <- est_Y1_a0$Onestep$estimated_psi
  onestep_Y1_a1_leq_EIF <- est_Y1_a1$Onestep$EIF
  onestep_Y1_a0_leq_EIF <- est_Y1_a0$Onestep$EIF
  
  tmle_Y1_a1 <- tmle_Y1_a1_leq - tmle_Y0_a1
  tmle_Y1_a0 <- tmle_Y1_a0_leq - tmle_Y0_a0
  tmle_Y1_a1_EIF <- tmle_Y1_a1_leq_EIF - tmle_Y0_a1_EIF
  tmle_Y1_a0_EIF <- tmle_Y1_a0_leq_EIF - tmle_Y0_a0_EIF
  onestep_Y1_a1 <- onestep_Y1_a1_leq - onestep_Y0_a1
  onestep_Y1_a0 <- onestep_Y1_a0_leq - onestep_Y0_a0
  onestep_Y1_a1_EIF <- onestep_Y1_a1_leq_EIF - onestep_Y0_a1_EIF
  onestep_Y1_a0_EIF <- onestep_Y1_a0_leq_EIF - onestep_Y0_a0_EIF
  
  ## get E[Y(a) <=2]
  #a=1
  est_Y2_a1 <- estfd(a=1,data=dt_Y2,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y2_a0 <- estfd(a=0,data=dt_Y2,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y2_a1_leq <- est_Y2_a1$TMLE$estimated_psi
  tmle_Y2_a0_leq <- est_Y2_a0$TMLE$estimated_psi
  tmle_Y2_a1_leq_EIF <- est_Y2_a1$TMLE$EIF
  tmle_Y2_a0_leq_EIF <- est_Y2_a0$TMLE$EIF
  onestep_Y2_a1_leq <- est_Y2_a1$Onestep$estimated_psi
  onestep_Y2_a0_leq <- est_Y2_a0$Onestep$estimated_psi
  onestep_Y2_a1_leq_EIF <- est_Y2_a1$Onestep$EIF
  onestep_Y2_a0_leq_EIF <- est_Y2_a0$Onestep$EIF
  
  tmle_Y2_a1 <- tmle_Y2_a1_leq - tmle_Y1_a1_leq
  tmle_Y2_a0 <- tmle_Y2_a0_leq - tmle_Y1_a0_leq
  tmle_Y2_a1_EIF <- tmle_Y2_a1_leq_EIF - tmle_Y1_a1_leq_EIF
  tmle_Y2_a0_EIF <- tmle_Y2_a0_leq_EIF - tmle_Y1_a0_leq_EIF
  onestep_Y2_a1 <- onestep_Y2_a1_leq - onestep_Y1_a1_leq
  onestep_Y2_a0 <- onestep_Y2_a0_leq - onestep_Y1_a0_leq
  onestep_Y2_a1_EIF <- onestep_Y2_a1_leq_EIF - onestep_Y1_a1_leq_EIF
  onestep_Y2_a0_EIF <- onestep_Y2_a0_leq_EIF - onestep_Y1_a0_leq_EIF
  
  ## get E[Y(a)<=3]
  #a=1
  est_Y3_a1 <- estfd(a=1,data=dt_Y3,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y3_a0 <- estfd(a=0,data=dt_Y3,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y3_a1_leq <- est_Y3_a1$TMLE$estimated_psi
  tmle_Y3_a0_leq <- est_Y3_a0$TMLE$estimated_psi
  tmle_Y3_a1_leq_EIF <- est_Y3_a1$TMLE$EIF
  tmle_Y3_a0_leq_EIF <- est_Y3_a0$TMLE$EIF
  onestep_Y3_a1_leq <- est_Y3_a1$Onestep$estimated_psi
  onestep_Y3_a0_leq <- est_Y3_a0$Onestep$estimated_psi
  onestep_Y3_a1_leq_EIF <- est_Y3_a1$Onestep$EIF
  onestep_Y3_a0_leq_EIF <- est_Y3_a0$Onestep$EIF
  
  tmle_Y3_a1 <- tmle_Y3_a1_leq - tmle_Y2_a1_leq
  tmle_Y3_a0 <- tmle_Y3_a0_leq - tmle_Y2_a0_leq
  tmle_Y3_a1_EIF <- tmle_Y3_a1_leq_EIF - tmle_Y2_a1_leq_EIF
  tmle_Y3_a0_EIF <- tmle_Y3_a0_leq_EIF - tmle_Y2_a0_leq_EIF
  onestep_Y3_a1 <- onestep_Y3_a1_leq - onestep_Y2_a1_leq
  onestep_Y3_a0 <- onestep_Y3_a0_leq - onestep_Y2_a0_leq
  onestep_Y3_a1_EIF <- onestep_Y3_a1_leq_EIF - onestep_Y2_a1_leq_EIF
  onestep_Y3_a0_EIF <- onestep_Y3_a0_leq_EIF - onestep_Y2_a0_leq_EIF
  
  ## get E[Y(a)<=4]
  #a=1
  est_Y4_a1 <- estfd(a=1,data=dt_Y4,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y4_a0 <- estfd(a=0,data=dt_Y4,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y4_a1_leq <- est_Y4_a1$TMLE$estimated_psi
  tmle_Y4_a0_leq <- est_Y4_a0$TMLE$estimated_psi
  tmle_Y4_a1_leq_EIF <- est_Y4_a1$TMLE$EIF
  tmle_Y4_a0_leq_EIF <- est_Y4_a0$TMLE$EIF
  onestep_Y4_a1_leq <- est_Y4_a1$Onestep$estimated_psi
  onestep_Y4_a0_leq <- est_Y4_a0$Onestep$estimated_psi
  onestep_Y4_a1_leq_EIF <- est_Y4_a1$Onestep$EIF
  onestep_Y4_a0_leq_EIF <- est_Y4_a0$Onestep$EIF
  
  tmle_Y4_a1 <- tmle_Y4_a1_leq - tmle_Y3_a1_leq
  tmle_Y4_a0 <- tmle_Y4_a0_leq - tmle_Y3_a0_leq
  tmle_Y4_a1_EIF <- tmle_Y4_a1_leq_EIF - tmle_Y3_a1_leq_EIF
  tmle_Y4_a0_EIF <- tmle_Y4_a0_leq_EIF - tmle_Y3_a0_leq_EIF
  onestep_Y4_a1 <- onestep_Y4_a1_leq - onestep_Y3_a1_leq
  onestep_Y4_a0 <- onestep_Y4_a0_leq - onestep_Y3_a0_leq
  onestep_Y4_a1_EIF <- onestep_Y4_a1_leq_EIF - onestep_Y3_a1_leq_EIF
  onestep_Y4_a0_EIF <- onestep_Y4_a0_leq_EIF - onestep_Y3_a0_leq_EIF
  
  ## get E[Y(a)<=5]
  #a=1
  est_Y5_a1 <- estfd(a=1,data=dt_Y5,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y5_a0 <- estfd(a=0,data=dt_Y5,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  tmle_Y5_a1_leq <- est_Y5_a1$TMLE$estimated_psi
  tmle_Y5_a0_leq <- est_Y5_a0$TMLE$estimated_psi
  tmle_Y5_a1_leq_EIF <- est_Y5_a1$TMLE$EIF
  tmle_Y5_a0_leq_EIF <- est_Y5_a0$TMLE$EIF
  onestep_Y5_a1_leq <- est_Y5_a1$Onestep$estimated_psi
  onestep_Y5_a0_leq <- est_Y5_a0$Onestep$estimated_psi
  onestep_Y5_a1_leq_EIF <- est_Y5_a1$Onestep$EIF
  onestep_Y5_a0_leq_EIF <- est_Y5_a0$Onestep$EIF
  
  tmle_Y5_a1 <- tmle_Y5_a1_leq - tmle_Y4_a1_leq
  tmle_Y5_a0 <- tmle_Y5_a0_leq - tmle_Y4_a0_leq
  tmle_Y5_a1_EIF <- tmle_Y5_a1_leq_EIF - tmle_Y4_a1_leq_EIF
  tmle_Y5_a0_EIF <- tmle_Y5_a0_leq_EIF - tmle_Y4_a0_leq_EIF
  onestep_Y5_a1 <- onestep_Y5_a1_leq - onestep_Y4_a1_leq
  onestep_Y5_a0 <- onestep_Y5_a0_leq - onestep_Y4_a0_leq
  onestep_Y5_a1_EIF <- onestep_Y5_a1_leq_EIF - onestep_Y4_a1_leq_EIF
  onestep_Y5_a0_EIF <- onestep_Y5_a0_leq_EIF - onestep_Y4_a0_leq_EIF
  
  ## get E[Y(a)==6]
  tmle_Y6_a1 <- 1 - tmle_Y5_a1 -tmle_Y4_a1 - tmle_Y3_a1 - tmle_Y2_a1 - tmle_Y1_a1 - tmle_Y0_a1
  tmle_Y6_a0 <- 1 - tmle_Y5_a0 -tmle_Y4_a0 - tmle_Y3_a0 - tmle_Y2_a0 - tmle_Y1_a0 - tmle_Y0_a0
  tmle_Y6_a1_EIF <- 1 - tmle_Y5_a1_EIF -tmle_Y4_a1_EIF - tmle_Y3_a1_EIF - tmle_Y2_a1_EIF - tmle_Y1_a1_EIF - tmle_Y0_a1_EIF
  tmle_Y6_a0_EIF <- 1 - tmle_Y5_a0_EIF -tmle_Y4_a0_EIF - tmle_Y3_a0_EIF - tmle_Y2_a0_EIF - tmle_Y1_a0_EIF - tmle_Y0_a0_EIF
  onestep_Y6_a1 <- 1 - onestep_Y5_a1 -onestep_Y4_a1 - onestep_Y3_a1 - onestep_Y2_a1 - onestep_Y1_a1 - onestep_Y0_a1
  onestep_Y6_a0 <- 1 - onestep_Y5_a0 -onestep_Y4_a0 - onestep_Y3_a0 - onestep_Y2_a0 - onestep_Y1_a0 - onestep_Y0_a0
  onestep_Y6_a1_EIF <- 1 - onestep_Y5_a1_EIF -onestep_Y4_a1_EIF - onestep_Y3_a1_EIF - onestep_Y2_a1_EIF - onestep_Y1_a1_EIF - onestep_Y0_a1_EIF
  onestep_Y6_a0_EIF <- 1 - onestep_Y5_a0_EIF -onestep_Y4_a0_EIF - onestep_Y3_a0_EIF - onestep_Y2_a0_EIF - onestep_Y1_a0_EIF - onestep_Y0_a0_EIF
  
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
  dim_tmle_a1_EIF <- 0*tmle_Y0_a1_EIF + 1*tmle_Y1_a1_EIF + 2*tmle_Y2_a1_EIF + 3*tmle_Y3_a1_EIF + 4*tmle_Y4_a1_EIF + 5*tmle_Y5_a1_EIF + 6*tmle_Y6_a1_EIF
  dim_tmle_a0_EIF <- 0*tmle_Y0_a0_EIF + 1*tmle_Y1_a0_EIF + 2*tmle_Y2_a0_EIF + 3*tmle_Y3_a0_EIF + 4*tmle_Y4_a0_EIF + 5*tmle_Y5_a0_EIF + 6*tmle_Y6_a0_EIF
  dim_tmle <- dim_tmle_a1 - dim_tmle_a0
  dim_tmle_EIF <- dim_tmle_a1_EIF - dim_tmle_a0_EIF
  
  dim_onestep_a1 <- 0*onestep_Y0_a1 + 1*onestep_Y1_a1 + 2*onestep_Y2_a1 + 3*onestep_Y3_a1 + 4*onestep_Y4_a1 + 5*onestep_Y5_a1 + 6*onestep_Y6_a1
  dim_onestep_a0 <- 0*onestep_Y0_a0 + 1*onestep_Y1_a0 + 2*onestep_Y2_a0 + 3*onestep_Y3_a0 + 4*onestep_Y4_a0 + 5*onestep_Y5_a0 + 6*onestep_Y6_a0
  dim_onestep_a1_EIF <- 0*onestep_Y0_a1_EIF + 1*onestep_Y1_a1_EIF + 2*onestep_Y2_a1_EIF + 3*onestep_Y3_a1_EIF + 4*onestep_Y4_a1_EIF + 5*onestep_Y5_a1_EIF + 6*onestep_Y6_a1_EIF
  dim_onestep_a0_EIF <- 0*onestep_Y0_a0_EIF + 1*onestep_Y1_a0_EIF + 2*onestep_Y2_a0_EIF + 3*onestep_Y3_a0_EIF + 4*onestep_Y4_a0_EIF + 5*onestep_Y5_a0_EIF + 6*onestep_Y6_a0_EIF
  dim_onestep <- dim_onestep_a1 - dim_onestep_a0
  dim_onestep_EIF <- dim_onestep_a1_EIF - dim_onestep_a0_EIF
  
  # PIIE
  piie_tmle_a1 <- dt$Y - dim_tmle_a1
  piie_tmle_a0 <- dt$Y - dim_tmle_a0
  piie_tmle_a1_EIF <- dt$Y - mean(dt$Y) - dim_tmle_a1_EIF
  piie_tmle_a0_EIF <- dt$Y - mean(dt$Y) - dim_tmle_a0_EIF
  
  piie_onestep_a1 <- dt$Y - dim_onestep_a1
  piie_onestep_a0 <- dt$Y - dim_onestep_a0
  piie_onestep_a1_EIF <- dt$Y - mean(dt$Y) - dim_onestep_a1_EIF
  piie_onestep_a0_EIF <- dt$Y - mean(dt$Y) - dim_onestep_a0_EIF
  
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
  
  return(list(
    # difference in mean
    dim_tmle = dim_tmle,
    dim_tmle_EIF = dim_tmle_EIF,
    dim_onestep = dim_onestep,
    dim_onestep_EIF = dim_onestep_EIF,
    # PIIE
    piie_tmle_a1 = piie_tmle_a1,
    piie_tmle_a0 = piie_tmle_a0,
    piie_tmle_a1_EIF = piie_tmle_a1_EIF,
    piie_tmle_a0_EIF = piie_tmle_a0_EIF,
    piie_onestep_a1 = piie_onestep_a1,
    piie_onestep_a0 = piie_onestep_a0,
    piie_onestep_a1_EIF = piie_onestep_a1_EIF,
    piie_onestep_a0_EIF = piie_onestep_a0_EIF,
    #LOR
    lor_tmle = lor_tmle,
    lor_onestep = lor_onestep,
    # Individual estimates
    tmle_Y0_a1 = tmle_Y0_a1, tmle_Y1_a1 = tmle_Y1_a1, tmle_Y2_a1 = tmle_Y2_a1, tmle_Y3_a1 = tmle_Y3_a1, tmle_Y4_a1 = tmle_Y4_a1,
    tmle_Y5_a1 = tmle_Y5_a1, tmle_Y6_a1 = tmle_Y6_a1,
    tmle_Y0_a0 = tmle_Y0_a0, tmle_Y1_a0 = tmle_Y1_a0, tmle_Y2_a0 = tmle_Y2_a0, tmle_Y3_a0 = tmle_Y3_a0, tmle_Y4_a0 = tmle_Y4_a0,
    tmle_Y5_a0 = tmle_Y5_a0, tmle_Y6_a0 = tmle_Y6_a0,
    onestep_Y0_a1 = onestep_Y0_a1, onestep_Y1_a1 = onestep_Y1_a1, onestep_Y2_a1 = onestep_Y2_a1, onestep_Y3_a1 = onestep_Y3_a1,
    onestep_Y4_a1 = onestep_Y4_a1, onestep_Y5_a1 = onestep_Y5_a1, onestep_Y6_a1 = onestep_Y6_a1,
    onestep_Y0_a0 = onestep_Y0_a0, onestep_Y1_a0 = onestep_Y1_a0, onestep_Y2_a0 = onestep_Y2_a0, onestep_Y3_a0 = onestep_Y3_a0,
    onestep_Y4_a0 = onestep_Y4_a0, onestep_Y5_a0 = onestep_Y5_a0, onestep_Y6_a0 = onestep_Y6_a0,
    # EIFs
    tmle_Y0_a1_EIF = tmle_Y0_a1_EIF, tmle_Y1_a1_EIF = tmle_Y1_a1_EIF, tmle_Y2_a1_EIF = tmle_Y2_a1_EIF, tmle_Y3_a1_EIF = tmle_Y3_a1_EIF,
    tmle_Y4_a1_EIF = tmle_Y4_a1_EIF, tmle_Y5_a1_EIF = tmle_Y5_a1_EIF, tmle_Y6_a1_EIF = tmle_Y6_a1_EIF,
    tmle_Y0_a0_EIF = tmle_Y0_a0_EIF, tmle_Y1_a0_EIF = tmle_Y1_a0_EIF, tmle_Y2_a0_EIF = tmle_Y2_a0_EIF, tmle_Y3_a0_EIF = tmle_Y3_a0_EIF,
    tmle_Y4_a0_EIF = tmle_Y4_a0_EIF, tmle_Y5_a0_EIF = tmle_Y5_a0_EIF, tmle_Y6_a0_EIF = tmle_Y6_a0_EIF,
    onestep_Y0_a1_EIF = onestep_Y0_a1_EIF, onestep_Y1_a1_EIF = onestep_Y1_a1_EIF, onestep_Y2_a1_EIF = onestep_Y2_a1_EIF,
    onestep_Y3_a1_EIF = onestep_Y3_a1_EIF, onestep_Y4_a1_EIF = onestep_Y4_a1_EIF, onestep_Y5_a1_EIF = onestep_Y5_a1_EIF,
    onestep_Y6_a1_EIF = onestep_Y6_a1_EIF, onestep_Y0_a0_EIF = onestep_Y0_a0_EIF, onestep_Y1_a0_EIF = onestep_Y1_a0_EIF,
    onestep_Y2_a0_EIF = onestep_Y2_a0_EIF, onestep_Y3_a0_EIF = onestep_Y3_a0_EIF, onestep_Y4_a0_EIF = onestep_Y4_a0_EIF,
    onestep_Y5_a0_EIF = onestep_Y5_a0_EIF, onestep_Y6_a0_EIF = onestep_Y6_a0_EIF
  ))
  
  
  
}

cat('Imputation', m, 'of 10\n')

dt.analysis <- dt %>% filter(.imp==m)
n <- nrow(dt.analysis)

# bootstrap for DIM and LOR
estimation.m <- TMLE_m(dt.analysis)

set1.TMLE.dim.est_m <- estimation.m$dim_tmle
set1.TMLE.piie.a1.est_m <- estimation.m$piie_tmle_a1
set1.TMLE.piie.a0.est_m <- estimation.m$piie_tmle_a0
set1.TMLE.lor.est_m <- estimation.m$lor_tmle

set1.onestep.dim.est_m <- estimation.m$dim_onestep
set1.onestep.piie.a1.est_m <- estimation.m$piie_onestep_a1
set1.onestep.piie.a0.est_m <- estimation.m$piie_onestep_a0
set1.onestep.lor.est_m <- estimation.m$lor_onestep

set1.TMLE.dim.var_m <- mean(estimation.m$dim_tmle_EIF^2)/n
set1.TMLE.piie.a1.var_m <- mean(estimation.m$piie_tmle_a1_EIF^2)/n
set1.TMLE.piie.a0.var_m <- mean(estimation.m$piie_tmle_a0_EIF^2)/n

set1.onestep.dim.var_m <-mean(estimation.m$dim_onestep_EIF^2)/n
set1.onestep.piie.a1.var_m <- mean(estimation.m$piie_onestep_a1_EIF^2)/n
set1.onestep.piie.a0.var_m <- mean(estimation.m$piie_onestep_a0_EIF^2)/n

TMLE.PMF.Ya1_m <- c(m, estimation.m$tmle_Y0_a1 ,estimation.m$tmle_Y1_a1, estimation.m$tmle_Y2_a1, estimation.m$tmle_Y3_a1, estimation.m$tmle_Y4_a1, estimation.m$tmle_Y5_a1, estimation.m$tmle_Y6_a1)
TMLE.PMF.Ya0_m <- c(m, estimation.m$tmle_Y0_a0 ,estimation.m$tmle_Y1_a0, estimation.m$tmle_Y2_a0, estimation.m$tmle_Y3_a0, estimation.m$tmle_Y4_a0, estimation.m$tmle_Y5_a0, estimation.m$tmle_Y6_a0)
onestep.PMF.Ya1_m <- c(m, estimation.m$onestep_Y0_a1 ,estimation.m$onestep_Y1_a1, estimation.m$onestep_Y2_a1, estimation.m$onestep_Y3_a1, estimation.m$onestep_Y4_a1, estimation.m$onestep_Y5_a1, estimation.m$onestep_Y6_a1)
onestep.PMF.Ya0_m <- c(m, estimation.m$onestep_Y0_a0 ,estimation.m$onestep_Y1_a0, estimation.m$onestep_Y2_a0, estimation.m$onestep_Y3_a0, estimation.m$onestep_Y4_a0, estimation.m$onestep_Y5_a0, estimation.m$onestep_Y6_a0)

save(list = c('estimation.m','set1.TMLE.dim.est_m','set1.TMLE.piie.a1.est_m','set1.TMLE.piie.a0.est_m',
              'set1.TMLE.lor.est_m','set1.onestep.dim.est_m','set1.onestep.piie.a1.est_m',
              'set1.onestep.piie.a0.est_m','set1.onestep.lor.est_m',
              'set1.TMLE.dim.var_m','set1.TMLE.piie.a1.var_m','set1.TMLE.piie.a0.var_m',
              'set1.onestep.dim.var_m','set1.onestep.piie.a1.var_m','set1.onestep.piie.a0.var_m',
              'TMLE.PMF.Ya1_m', 'TMLE.PMF.Ya0_m', 'onestep.PMF.Ya1_m', 'onestep.PMF.Ya0_m'),
     file = paste0('output/output_', m, '.RData'))

