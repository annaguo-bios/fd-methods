args = commandArgs(trailingOnly=T)
m=as.integer(args[1]) # m-th multiple imputation
set=as.integer(args[2]) # set of analysis (1,2,3,4,5,6,7)

set.seed(7)

library(fdcausal)
library(SuperLearner)
library(dplyr)
library(mice)
library(boot)


if(set==1){
  m.var=c('M1','M2')
  y.var='Y'
}else if(set==2){
  m.var="M.cutoff1"
  y.var='Y.cutoff2'
}else if(set==3){
  m.var="M.cutoff1"
  y.var='Y.cutoff3'
}else if(set==4){
  m.var="M.cutoff1"
  y.var='Y.cutoff4'
}else if(set==5){
  m.var="M.cutoff2"
  y.var='Y.cutoff2'
}else if(set==6){
  m.var="M.cutoff2"
  y.var='Y.cutoff3'
}else if(set==7){
  m.var="M.cutoff2"
  y.var='Y.cutoff4'
}

if(set==1){
  
  dt2 <- read.csv("export_continuousM.csv") %>% mutate(mds_alarm_to_needle = if_else(mds_tpa == 0, 0, mds_alarm_to_needle))
  dt2.imputed <- mice(dt2, m = 10, method = 'pmm', seed = 7)
  dt2.imputed_data <- complete(dt2.imputed, 'long')
  
  dt <- dt2.imputed_data %>% 
    # select(-X) %>% 
    rename(A = STEMO, Y = mrs, M = mds_alarm_to_needle, X1 = systolic_before, X2 = nihss_before) %>% 
    mutate(A = as.numeric(A), Y = as.numeric(Y)) %>% 
    mutate(M1 = mds_tpa, M2=as.numeric(mds_tpa*M))
  
}else{
  
  dt <- read.csv("export_mi.csv") %>% 
    select(-X) %>% 
    rename(A = STEMO, Y = mrs, M = mds_alarm_to_needle_cat, X1 = systolic_before, X2 = nihss_before) %>% 
    mutate(A = as.numeric(A), Y = as.numeric(Y)) %>% 
    mutate(M = case_when(
      M == '[20,48]' ~ 1,
      M == '(48,75]' ~ 2,
      M == '(75,Inf]' ~ 3
    )) %>%
    mutate(M = as.numeric(M)) %>%
    mutate(M.cutoff1 = 1 - as.numeric(M<=1), M.cutoff2 = 1 - as.numeric(M<=2)) %>% # binarize mediator M
    mutate(Y.cutoff2 = 1- as.numeric(Y<=2), Y.cutoff3 = 1 - as.numeric(Y<=3), Y.cutoff4 = 1 - as.numeric(Y<=4)) # binarize outcome Y
  
  
}



cat("Imputation", m, "\n")

dt.analysis <- dt %>% filter(.imp==m)
n <- nrow(dt.analysis)

est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=m.var,
            outcome=y.var, covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
            lib = c("SL.glm", "SL.earth" ,"SL.ranger",  "SL.mean"), cvg.criteria = n^{-0.5})

set.TMLE.est_m <- est$TMLE$ATE
set.onestep.est_m <- est$Onestep$ATE
set.piie.a1.TMLE.est_m <- mean(dt.analysis[[y.var]]) - est$TMLE.Y1$estimated_psi
set.piie.a0.TMLE.est_m <- mean(dt.analysis[[y.var]]) - est$TMLE.Y0$estimated_psi
set.piie.a1.onestep.est_m <- mean(dt.analysis[[y.var]]) - est$Onestep.Y1$estimated_psi
set.piie.a0.onestep.est_m <- mean(dt.analysis[[y.var]]) - est$Onestep.Y0$estimated_psi

set.TMLE.var_m <- mean(est$TMLE$EIF^2)/n
set.onestep.var_m <- mean(est$Onestep$EIF^2)/n
set.piie.a1.TMLE.var_m <- mean({(dt.analysis[[y.var]] - mean(dt.analysis[[y.var]])) - est$TMLE.Y1$EIF}^2)/n
set.piie.a0.TMLE.var_m <- mean({(dt.analysis[[y.var]] - mean(dt.analysis[[y.var]])) - est$TMLE.Y0$EIF}^2)/n
set.piie.a1.onestep.var_m <- mean({(dt.analysis[[y.var]] - mean(dt.analysis[[y.var]])) - est$Onestep.Y1$EIF}^2)/n
set.piie.a0.onestep.var_m <- mean({(dt.analysis[[y.var]] - mean(dt.analysis[[y.var]])) - est$Onestep.Y0$EIF}^2)/n

save(list=c("est","set.TMLE.est_m", "set.onestep.est_m", "set.piie.a1.TMLE.est_m", "set.piie.a0.TMLE.est_m",
     "set.piie.a1.onestep.est_m", "set.piie.a0.onestep.est_m",
     "set.TMLE.var_m", "set.onestep.var_m", "set.piie.a1.TMLE.var_m", "set.piie.a0.TMLE.var_m",
     "set.piie.a1.onestep.var_m", "set.piie.a0.onestep.var_m"),
     file=paste0("output/set_", set, "_m_", m, ".RData"))
