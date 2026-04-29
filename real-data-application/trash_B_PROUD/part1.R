library(fdcausal)
library(SuperLearner)
library(dplyr)
library(mice)
library(boot)

# Binarize M & Y or treat them as continuous ========================

## compare whether the original data (dt2) with continuous M agree with the organized categorized data (dt1)
dt1 <- read.csv("../../data/B_PROUD/export_mi.csv") %>% filter(.imp==0)
dt2 <- read.csv("../../data/B_PROUD/export_continuousM.csv") %>% mutate(mds_alarm_to_needle = if_else(mds_tpa == 0, 0, mds_alarm_to_needle))
dt2.imputed <- mice(dt2, m = 10, method = 'pmm', seed = 7)
dt2.imputed_data <- complete(dt2.imputed, 'long')

dt <- dt2.imputed_data %>% 
  # select(-X) %>% 
  rename(A = STEMO, Y = mrs, M = mds_alarm_to_needle, X1 = systolic_before, X2 = nihss_before) %>% 
  mutate(A = as.numeric(A), Y = as.numeric(Y)) %>% 
  mutate(M1 = mds_tpa, M2=as.numeric(mds_tpa*M))



## Setting 1: M, Y as continuous ====
set1.TMLE.est <- rep(NA,10); set1.TMLE.var <- rep(NA,10)
set1.onestep.est <- rep(NA,10); set1.onestep.var <- rep(NA,10)

set1.piie.a1.TMLE.est <- rep(NA,10); set1.piie.a1.TMLE.var <- rep(NA,10)
set1.piie.a0.TMLE.est <- rep(NA,10); set1.piie.a0.TMLE.var <- rep(NA,10)
set1.piie.a1.onestep.est <- rep(NA,10); set1.piie.a1.onestep.var <- rep(NA,10)
set1.piie.a0.onestep.est <- rep(NA,10); set1.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Imputation", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M1","M2"),
              outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set1.TMLE.est[m] <- est$TMLE$ATE
  set1.onestep.est[m] <- est$Onestep$ATE
  set1.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y) - est$TMLE.Y1$estimated_psi
  set1.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y) - est$TMLE.Y0$estimated_psi
  set1.piie.a1.onestep.est[m] <- mean(dt.analysis$Y) - est$Onestep.Y1$estimated_psi
  set1.piie.a0.onestep.est[m] <- mean(dt.analysis$Y) - est$Onestep.Y0$estimated_psi
  
  set1.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set1.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set1.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y - mean(dt.analysis$Y)) - est$TMLE.Y1$EIF}^2)/n
  set1.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y - mean(dt.analysis$Y)) - est$TMLE.Y0$EIF}^2)/n
  set1.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y - mean(dt.analysis$Y)) - est$Onestep.Y1$EIF}^2)/n
  set1.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y - mean(dt.analysis$Y)) - est$Onestep.Y0$EIF}^2)/n
}

set1.TMLE.mean <- mean(set1.TMLE.est)
set1.onestep.mean <- mean(set1.onestep.est)
set1.piie.a1.TMLE.mean <- mean(set1.piie.a1.TMLE.est)
set1.piie.a0.TMLE.mean <- mean(set1.piie.a0.TMLE.est)
set1.piie.a1.onestep.mean <- mean(set1.piie.a1.onestep.est)
set1.piie.a0.onestep.mean <- mean(set1.piie.a0.onestep.est)

set1.TMLE.imp.var <- mean(set1.TMLE.var) + sum((set1.TMLE.est - set1.TMLE.mean)^2)/9
set1.onestep.imp.var <- mean(set1.onestep.var) + sum((set1.onestep.est - set1.onestep.mean)^2)/9
set1.piie.a1.TMLE.imp.var <- mean(set1.piie.a1.TMLE.var) + sum((set1.piie.a1.TMLE.est - set1.piie.a1.TMLE.mean)^2)/9
set1.piie.a0.TMLE.imp.var <- mean(set1.piie.a0.TMLE.var) + sum((set1.piie.a0.TMLE.est - set1.piie.a0.TMLE.mean)^2)/9
set1.piie.a1.onestep.imp.var <- mean(set1.piie.a1.onestep.var) + sum((set1.piie.a1.onestep.est - set1.piie.a1.onestep.mean)^2)/9
set1.piie.a0.onestep.imp.var <- mean(set1.piie.a0.onestep.var) + sum((set1.piie.a0.onestep.est - set1.piie.a0.onestep.mean)^2)/9

round(set1.TMLE.mean - 1.96*sqrt(set1.TMLE.imp.var),3); round(set1.TMLE.mean + 1.96*sqrt(set1.TMLE.imp.var),3)
round(set1.onestep.mean - 1.96*sqrt(set1.onestep.imp.var),3); round(set1.onestep.mean + 1.96*sqrt(set1.onestep.imp.var),3)
round(set1.piie.a1.TMLE.mean - 1.96*sqrt(set1.piie.a1.TMLE.imp.var),3); round(set1.piie.a1.TMLE.mean + 1.96*sqrt(set1.piie.a1.TMLE.imp.var),3)
round(set1.piie.a0.TMLE.mean - 1.96*sqrt(set1.piie.a0.TMLE.imp.var),3); round(set1.piie.a0.TMLE.mean + 1.96*sqrt(set1.piie.a0.TMLE.imp.var),3)
round(set1.piie.a1.onestep.mean - 1.96*sqrt(set1.piie.a1.onestep.imp.var),3); round(set1.piie.a1.onestep.mean + 1.96*sqrt(set1.piie.a1.onestep.imp.var),3)
round(set1.piie.a0.onestep.mean - 1.96*sqrt(set1.piie.a0.onestep.imp.var),3); round(set1.piie.a0.onestep.mean + 1.96*sqrt(set1.piie.a0.onestep.imp.var),3)

set1.result <- list(TMLE = data.frame(est = set1.TMLE.mean, var= set1.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set1.piie.a1.TMLE.mean, var= set1.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set1.piie.a0.TMLE.mean, var= set1.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set1.onestep.mean, var=set1.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set1.piie.a1.onestep.mean, var= set1.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set1.piie.a0.onestep.mean, var= set1.piie.a0.onestep.imp.var)
)

## Setting 2: M.cutoff1, Y.cutoff2 ====
dt <- read.csv("../../data/B_PROUD/export_mi.csv") %>% 
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

set2.TMLE.est <- rep(NA,10); set2.TMLE.var <- rep(NA,10)
set2.onestep.est <- rep(NA,10); set2.onestep.var <- rep(NA,10)
set2.piie.a1.TMLE.est <- rep(NA,10); set2.piie.a1.TMLE.var <- rep(NA,10)
set2.piie.a0.TMLE.est <- rep(NA,10); set2.piie.a0.TMLE.var <- rep(NA,10)
set2.piie.a1.onestep.est <- rep(NA,10); set2.piie.a1.onestep.var <- rep(NA,10)
set2.piie.a0.onestep.est <- rep(NA,10); set2.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Setting 2, Imp: ", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff1"),
              outcome="Y.cutoff2", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set2.TMLE.est[m] <- est$TMLE$ATE
  set2.onestep.est[m] <- est$Onestep$ATE
  set2.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff2) - est$TMLE.Y1$estimated_psi
  set2.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff2) - est$TMLE.Y0$estimated_psi
  set2.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff2) - est$Onestep.Y1$estimated_psi
  set2.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff2) - est$Onestep.Y0$estimated_psi
  
  set2.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set2.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set2.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$TMLE.Y1$EIF}^2)/n
  set2.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$TMLE.Y0$EIF}^2)/n
  set2.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$Onestep.Y1$EIF}^2)/n
  set2.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$Onestep.Y0$EIF}^2)/n
}

set2.TMLE.mean <- mean(set2.TMLE.est)
set2.onestep.mean <- mean(set2.onestep.est)
set2.piie.a1.TMLE.mean <- mean(set2.piie.a1.TMLE.est)
set2.piie.a0.TMLE.mean <- mean(set2.piie.a0.TMLE.est)
set2.piie.a1.onestep.mean <- mean(set2.piie.a1.onestep.est)
set2.piie.a0.onestep.mean <- mean(set2.piie.a0.onestep.est)

set2.TMLE.imp.var <- mean(set2.TMLE.var) + sum((set2.TMLE.est - set2.TMLE.mean)^2)/9
set2.onestep.imp.var <- mean(set2.onestep.var) + sum((set2.onestep.est - set2.onestep.mean)^2)/9
set2.piie.a1.TMLE.imp.var <- mean(set2.piie.a1.TMLE.var) + sum((set2.piie.a1.TMLE.est - set2.piie.a1.TMLE.mean)^2)/9
set2.piie.a0.TMLE.imp.var <- mean(set2.piie.a0.TMLE.var) + sum((set2.piie.a0.TMLE.est - set2.piie.a0.TMLE.mean)^2)/9
set2.piie.a1.onestep.imp.var <- mean(set2.piie.a1.onestep.var) + sum((set2.piie.a1.onestep.est - set2.piie.a1.onestep.mean)^2)/9
set2.piie.a0.onestep.imp.var <- mean(set2.piie.a0.onestep.var) + sum((set2.piie.a0.onestep.est - set2.piie.a0.onestep.mean)^2)/9

round(set2.TMLE.mean - 1.96*sqrt(set2.TMLE.imp.var),3); round(set2.TMLE.mean + 1.96*sqrt(set2.TMLE.imp.var),3)
round(set2.onestep.mean - 1.96*sqrt(set2.onestep.imp.var),3); round(set2.onestep.mean + 1.96*sqrt(set2.onestep.imp.var),3)
round(set2.piie.a1.TMLE.mean - 1.96*sqrt(set2.piie.a1.TMLE.imp.var),3); round(set2.piie.a1.TMLE.mean + 1.96*sqrt(set2.piie.a1.TMLE.imp.var),3)
round(set2.piie.a0.TMLE.mean - 1.96*sqrt(set2.piie.a0.TMLE.imp.var),3); round(set2.piie.a0.TMLE.mean + 1.96*sqrt(set2.piie.a0.TMLE.imp.var),3)
round(set2.piie.a1.onestep.mean - 1.96*sqrt(set2.piie.a1.onestep.imp.var),3); round(set2.piie.a1.onestep.mean + 1.96*sqrt(set2.piie.a1.onestep.imp.var),3)
round(set2.piie.a0.onestep.mean - 1.96*sqrt(set2.piie.a0.onestep.imp.var),3); round(set2.piie.a0.onestep.mean + 1.96*sqrt(set2.piie.a0.onestep.imp.var),3)

set2.result <- list(TMLE = data.frame(est = set2.TMLE.mean, var= set2.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set2.piie.a1.TMLE.mean, var= set2.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set2.piie.a0.TMLE.mean, var= set2.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set2.onestep.mean, var=set2.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set2.piie.a1.onestep.mean, var= set2.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set2.piie.a0.onestep.mean, var= set2.piie.a0.onestep.imp.var)
)


## Setting 3: M.cutoff1, Y.cutoff3 ====
set3.TMLE.est <- rep(NA,10); set3.TMLE.var <- rep(NA,10)
set3.onestep.est <- rep(NA,10); set3.onestep.var <- rep(NA,10)
set3.piie.a1.TMLE.est <- rep(NA,10); set3.piie.a1.TMLE.var <- rep(NA,10)
set3.piie.a0.TMLE.est <- rep(NA,10); set3.piie.a0.TMLE.var <- rep(NA,10)
set3.piie.a1.onestep.est <- rep(NA,10); set3.piie.a1.onestep.var <- rep(NA,10)
set3.piie.a0.onestep.est <- rep(NA,10); set3.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Setting 3, Imp: ", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff1"),
              outcome="Y.cutoff3", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set3.TMLE.est[m] <- est$TMLE$ATE
  set3.onestep.est[m] <- est$Onestep$ATE
  set3.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff3) - est$TMLE.Y1$estimated_psi
  set3.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff3) - est$TMLE.Y0$estimated_psi
  set3.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff3) - est$Onestep.Y1$estimated_psi
  set3.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff3) - est$Onestep.Y0$estimated_psi
  
  set3.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set3.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set3.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$TMLE.Y1$EIF}^2)/n
  set3.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$TMLE.Y0$EIF}^2)/n
  set3.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$Onestep.Y1$EIF}^2)/n
  set3.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$Onestep.Y0$EIF}^2)/n
}

set3.TMLE.mean <- mean(set3.TMLE.est)
set3.onestep.mean <- mean(set3.onestep.est)
set3.piie.a1.TMLE.mean <- mean(set3.piie.a1.TMLE.est)
set3.piie.a0.TMLE.mean <- mean(set3.piie.a0.TMLE.est)
set3.piie.a1.onestep.mean <- mean(set3.piie.a1.onestep.est)
set3.piie.a0.onestep.mean <- mean(set3.piie.a0.onestep.est)

set3.TMLE.imp.var <- mean(set3.TMLE.var) + sum((set3.TMLE.est - set3.TMLE.mean)^2)/9
set3.onestep.imp.var <- mean(set3.onestep.var) + sum((set3.onestep.est - set3.onestep.mean)^2)/9
set3.piie.a1.TMLE.imp.var <- mean(set3.piie.a1.TMLE.var) + sum((set3.piie.a1.TMLE.est - set3.piie.a1.TMLE.mean)^2)/9
set3.piie.a0.TMLE.imp.var <- mean(set3.piie.a0.TMLE.var) + sum((set3.piie.a0.TMLE.est - set3.piie.a0.TMLE.mean)^2)/9
set3.piie.a1.onestep.imp.var <- mean(set3.piie.a1.onestep.var) + sum((set3.piie.a1.onestep.est - set3.piie.a1.onestep.mean)^2)/9
set3.piie.a0.onestep.imp.var <- mean(set3.piie.a0.onestep.var) + sum((set3.piie.a0.onestep.est - set3.piie.a0.onestep.mean)^2)/9

round(set3.TMLE.mean - 1.96*sqrt(set3.TMLE.imp.var),3); round(set3.TMLE.mean + 1.96*sqrt(set3.TMLE.imp.var),3)
round(set3.onestep.mean - 1.96*sqrt(set3.onestep.imp.var),3); round(set3.onestep.mean + 1.96*sqrt(set3.onestep.imp.var),3)
round(set3.piie.a1.TMLE.mean - 1.96*sqrt(set3.piie.a1.TMLE.imp.var),3); round(set3.piie.a1.TMLE.mean + 1.96*sqrt(set3.piie.a1.TMLE.imp.var),3)
round(set3.piie.a0.TMLE.mean - 1.96*sqrt(set3.piie.a0.TMLE.imp.var),3); round(set3.piie.a0.TMLE.mean + 1.96*sqrt(set3.piie.a0.TMLE.imp.var),3)
round(set3.piie.a1.onestep.mean - 1.96*sqrt(set3.piie.a1.onestep.imp.var),3); round(set3.piie.a1.onestep.mean + 1.96*sqrt(set3.piie.a1.onestep.imp.var),3)
round(set3.piie.a0.onestep.mean - 1.96*sqrt(set3.piie.a0.onestep.imp.var),3); round(set3.piie.a0.onestep.mean + 1.96*sqrt(set3.piie.a0.onestep.imp.var),3)

set3.result <- list(TMLE = data.frame(est = set3.TMLE.mean, var= set3.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set3.piie.a1.TMLE.mean, var= set3.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set3.piie.a0.TMLE.mean, var= set3.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set3.onestep.mean, var=set3.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set3.piie.a1.onestep.mean, var= set3.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set3.piie.a0.onestep.mean, var= set3.piie.a0.onestep.imp.var)
)

## Setting 4: M.cutoff1, Y.cutoff4 ====
set4.TMLE.est <- rep(NA,10); set4.TMLE.var <- rep(NA,10)
set4.onestep.est <- rep(NA,10); set4.onestep.var <- rep(NA,10)
set4.piie.a1.TMLE.est <- rep(NA,10); set4.piie.a1.TMLE.var <- rep(NA,10)
set4.piie.a0.TMLE.est <- rep(NA,10); set4.piie.a0.TMLE.var <- rep(NA,10)
set4.piie.a1.onestep.est <- rep(NA,10); set4.piie.a1.onestep.var <- rep(NA,10)
set4.piie.a0.onestep.est <- rep(NA,10); set4.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Setting 4, Imp: ", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff1"),
              outcome="Y.cutoff4", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set4.TMLE.est[m] <- est$TMLE$ATE
  set4.onestep.est[m] <- est$Onestep$ATE
  set4.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff4) - est$TMLE.Y1$estimated_psi
  set4.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff4) - est$TMLE.Y0$estimated_psi
  set4.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff4) - est$Onestep.Y1$estimated_psi
  set4.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff4) - est$Onestep.Y0$estimated_psi
  
  set4.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set4.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set4.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$TMLE.Y1$EIF}^2)/n
  set4.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$TMLE.Y0$EIF}^2)/n
  set4.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$Onestep.Y1$EIF}^2)/n
  set4.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$Onestep.Y0$EIF}^2)/n
}

set4.TMLE.mean <- mean(set4.TMLE.est)
set4.onestep.mean <- mean(set4.onestep.est)
set4.piie.a1.TMLE.mean <- mean(set4.piie.a1.TMLE.est)
set4.piie.a0.TMLE.mean <- mean(set4.piie.a0.TMLE.est)
set4.piie.a1.onestep.mean <- mean(set4.piie.a1.onestep.est)
set4.piie.a0.onestep.mean <- mean(set4.piie.a0.onestep.est)

set4.TMLE.imp.var <- mean(set4.TMLE.var) + sum((set4.TMLE.est - set4.TMLE.mean)^2)/9
set4.onestep.imp.var <- mean(set4.onestep.var) + sum((set4.onestep.est - set4.onestep.mean)^2)/9
set4.piie.a1.TMLE.imp.var <- mean(set4.piie.a1.TMLE.var) + sum((set4.piie.a1.TMLE.est - set4.piie.a1.TMLE.mean)^2)/9
set4.piie.a0.TMLE.imp.var <- mean(set4.piie.a0.TMLE.var) + sum((set4.piie.a0.TMLE.est - set4.piie.a0.TMLE.mean)^2)/9
set4.piie.a1.onestep.imp.var <- mean(set4.piie.a1.onestep.var) + sum((set4.piie.a1.onestep.est - set4.piie.a1.onestep.mean)^2)/9
set4.piie.a0.onestep.imp.var <- mean(set4.piie.a0.onestep.var) + sum((set4.piie.a0.onestep.est - set4.piie.a0.onestep.mean)^2)/9

round(set4.TMLE.mean - 1.96*sqrt(set4.TMLE.imp.var),3); round(set4.TMLE.mean + 1.96*sqrt(set4.TMLE.imp.var),3)
round(set4.onestep.mean - 1.96*sqrt(set4.onestep.imp.var),3); round(set4.onestep.mean + 1.96*sqrt(set4.onestep.imp.var),3)
round(set4.piie.a1.TMLE.mean - 1.96*sqrt(set4.piie.a1.TMLE.imp.var),3); round(set4.piie.a1.TMLE.mean + 1.96*sqrt(set4.piie.a1.TMLE.imp.var),3)
round(set4.piie.a0.TMLE.mean - 1.96*sqrt(set4.piie.a0.TMLE.imp.var),3); round(set4.piie.a0.TMLE.mean + 1.96*sqrt(set4.piie.a0.TMLE.imp.var),3)
round(set4.piie.a1.onestep.mean - 1.96*sqrt(set4.piie.a1.onestep.imp.var),3); round(set4.piie.a1.onestep.mean + 1.96*sqrt(set4.piie.a1.onestep.imp.var),3)
round(set4.piie.a0.onestep.mean - 1.96*sqrt(set4.piie.a0.onestep.imp.var),3); round(set4.piie.a0.onestep.mean + 1.96*sqrt(set4.piie.a0.onestep.imp.var),3)

set4.result <- list(TMLE = data.frame(est = set4.TMLE.mean, var= set4.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set4.piie.a1.TMLE.mean, var= set4.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set4.piie.a0.TMLE.mean, var= set4.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set4.onestep.mean, var=set4.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set4.piie.a1.onestep.mean, var= set4.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set4.piie.a0.onestep.mean, var= set4.piie.a0.onestep.imp.var)
)

## Setting 5: M.cutoff2, Y.cutoff2 ====
set5.TMLE.est <- rep(NA,10); set5.TMLE.var <- rep(NA,10)
set5.onestep.est <- rep(NA,10); set5.onestep.var <- rep(NA,10)
set5.piie.a1.TMLE.est <- rep(NA,10); set5.piie.a1.TMLE.var <- rep(NA,10)
set5.piie.a0.TMLE.est <- rep(NA,10); set5.piie.a0.TMLE.var <- rep(NA,10)
set5.piie.a1.onestep.est <- rep(NA,10); set5.piie.a1.onestep.var <- rep(NA,10)
set5.piie.a0.onestep.est <- rep(NA,10); set5.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat(sprintf("Setting 5, Imp: %d", m))
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff2"),
              outcome="Y.cutoff2", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set5.TMLE.est[m] <- est$TMLE$ATE
  set5.onestep.est[m] <- est$Onestep$ATE
  set5.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff2) - est$TMLE.Y1$estimated_psi
  set5.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff2) - est$TMLE.Y0$estimated_psi
  set5.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff2) - est$Onestep.Y1$estimated_psi
  set5.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff2) - est$Onestep.Y0$estimated_psi
  
  set5.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set5.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set5.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$TMLE.Y1$EIF}^2)/n
  set5.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$TMLE.Y0$EIF}^2)/n
  set5.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$Onestep.Y1$EIF}^2)/n
  set5.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff2 - mean(dt.analysis$Y.cutoff2)) - est$Onestep.Y0$EIF}^2)/n
}

set5.TMLE.mean <- mean(set5.TMLE.est)
set5.onestep.mean <- mean(set5.onestep.est)
set5.piie.a1.TMLE.mean <- mean(set5.piie.a1.TMLE.est)
set5.piie.a0.TMLE.mean <- mean(set5.piie.a0.TMLE.est)
set5.piie.a1.onestep.mean <- mean(set5.piie.a1.onestep.est)
set5.piie.a0.onestep.mean <- mean(set5.piie.a0.onestep.est)

set5.TMLE.imp.var <- mean(set5.TMLE.var) + sum((set5.TMLE.est - set5.TMLE.mean)^2)/9
set5.onestep.imp.var <- mean(set5.onestep.var) + sum((set5.onestep.est - set5.onestep.mean)^2)/9
set5.piie.a1.TMLE.imp.var <- mean(set5.piie.a1.TMLE.var) + sum((set5.piie.a1.TMLE.est - set5.piie.a1.TMLE.mean)^2)/9
set5.piie.a0.TMLE.imp.var <- mean(set5.piie.a0.TMLE.var) + sum((set5.piie.a0.TMLE.est - set5.piie.a0.TMLE.mean)^2)/9
set5.piie.a1.onestep.imp.var <- mean(set5.piie.a1.onestep.var) + sum((set5.piie.a1.onestep.est - set5.piie.a1.onestep.mean)^2)/9
set5.piie.a0.onestep.imp.var <- mean(set5.piie.a0.onestep.var) + sum((set5.piie.a0.onestep.est - set5.piie.a0.onestep.mean)^2)/9

round(set5.TMLE.mean - 1.96*sqrt(set5.TMLE.imp.var),3); round(set5.TMLE.mean + 1.96*sqrt(set5.TMLE.imp.var),3)
round(set5.onestep.mean - 1.96*sqrt(set5.onestep.imp.var),3); round(set5.onestep.mean + 1.96*sqrt(set5.onestep.imp.var),3)
round(set5.piie.a1.TMLE.mean - 1.96*sqrt(set5.piie.a1.TMLE.imp.var),3); round(set5.piie.a1.TMLE.mean + 1.96*sqrt(set5.piie.a1.TMLE.imp.var),3)
round(set5.piie.a0.TMLE.mean - 1.96*sqrt(set5.piie.a0.TMLE.imp.var),3); round(set5.piie.a0.TMLE.mean + 1.96*sqrt(set5.piie.a0.TMLE.imp.var),3)
round(set5.piie.a1.onestep.mean - 1.96*sqrt(set5.piie.a1.onestep.imp.var),3); round(set5.piie.a1.onestep.mean + 1.96*sqrt(set5.piie.a1.onestep.imp.var),3)
round(set5.piie.a0.onestep.mean - 1.96*sqrt(set5.piie.a0.onestep.imp.var),3); round(set5.piie.a0.onestep.mean + 1.96*sqrt(set5.piie.a0.onestep.imp.var),3)

set5.result <- list(TMLE = data.frame(est = set5.TMLE.mean, var= set5.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set5.piie.a1.TMLE.mean, var= set5.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set5.piie.a0.TMLE.mean, var= set5.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set5.onestep.mean, var=set5.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set5.piie.a1.onestep.mean, var= set5.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set5.piie.a0.onestep.mean, var= set5.piie.a0.onestep.imp.var)
)

## Setting 6: M.cutoff2, Y.cutoff3 ====
set6.TMLE.est <- rep(NA,10); set6.TMLE.var <- rep(NA,10)
set6.onestep.est <- rep(NA,10); set6.onestep.var <- rep(NA,10)
set6.piie.a1.TMLE.est <- rep(NA,10); set6.piie.a1.TMLE.var <- rep(NA,10)
set6.piie.a0.TMLE.est <- rep(NA,10); set6.piie.a0.TMLE.var <- rep(NA,10)
set6.piie.a1.onestep.est <- rep(NA,10); set6.piie.a1.onestep.var <- rep(NA,10)
set6.piie.a0.onestep.est <- rep(NA,10); set6.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Setting 6, Imp: ", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff2"),
              outcome="Y.cutoff3", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set6.TMLE.est[m] <- est$TMLE$ATE
  set6.onestep.est[m] <- est$Onestep$ATE
  set6.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff3) - est$TMLE.Y1$estimated_psi
  set6.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff3) - est$TMLE.Y0$estimated_psi
  set6.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff3) - est$Onestep.Y1$estimated_psi
  set6.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff3) - est$Onestep.Y0$estimated_psi
  
  set6.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set6.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set6.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$TMLE.Y1$EIF}^2)/n
  set6.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$TMLE.Y0$EIF}^2)/n
  set6.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$Onestep.Y1$EIF}^2)/n
  set6.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff3 - mean(dt.analysis$Y.cutoff3)) - est$Onestep.Y0$EIF}^2)/n
}

set6.TMLE.mean <- mean(set6.TMLE.est)
set6.onestep.mean <- mean(set6.onestep.est)
set6.piie.a1.TMLE.mean <- mean(set6.piie.a1.TMLE.est)
set6.piie.a0.TMLE.mean <- mean(set6.piie.a0.TMLE.est)
set6.piie.a1.onestep.mean <- mean(set6.piie.a1.onestep.est)
set6.piie.a0.onestep.mean <- mean(set6.piie.a0.onestep.est)

set6.TMLE.imp.var <- mean(set6.TMLE.var) + sum((set6.TMLE.est - set6.TMLE.mean)^2)/9
set6.onestep.imp.var <- mean(set6.onestep.var) + sum((set6.onestep.est - set6.onestep.mean)^2)/9
set6.piie.a1.TMLE.imp.var <- mean(set6.piie.a1.TMLE.var) + sum((set6.piie.a1.TMLE.est - set6.piie.a1.TMLE.mean)^2)/9
set6.piie.a0.TMLE.imp.var <- mean(set6.piie.a0.TMLE.var) + sum((set6.piie.a0.TMLE.est - set6.piie.a0.TMLE.mean)^2)/9
set6.piie.a1.onestep.imp.var <- mean(set6.piie.a1.onestep.var) + sum((set6.piie.a1.onestep.est - set6.piie.a1.onestep.mean)^2)/9
set6.piie.a0.onestep.imp.var <- mean(set6.piie.a0.onestep.var) + sum((set6.piie.a0.onestep.est - set6.piie.a0.onestep.mean)^2)/9

round(set6.TMLE.mean - 1.96*sqrt(set6.TMLE.imp.var),3); round(set6.TMLE.mean + 1.96*sqrt(set6.TMLE.imp.var),3)
round(set6.onestep.mean - 1.96*sqrt(set6.onestep.imp.var),3); round(set6.onestep.mean + 1.96*sqrt(set6.onestep.imp.var),3)
round(set6.piie.a1.TMLE.mean - 1.96*sqrt(set6.piie.a1.TMLE.imp.var),3); round(set6.piie.a1.TMLE.mean + 1.96*sqrt(set6.piie.a1.TMLE.imp.var),3)
round(set6.piie.a0.TMLE.mean - 1.96*sqrt(set6.piie.a0.TMLE.imp.var),3); round(set6.piie.a0.TMLE.mean + 1.96*sqrt(set6.piie.a0.TMLE.imp.var),3)
round(set6.piie.a1.onestep.mean - 1.96*sqrt(set6.piie.a1.onestep.imp.var),3); round(set6.piie.a1.onestep.mean + 1.96*sqrt(set6.piie.a1.onestep.imp.var),3)
round(set6.piie.a0.onestep.mean - 1.96*sqrt(set6.piie.a0.onestep.imp.var),3); round(set6.piie.a0.onestep.mean + 1.96*sqrt(set6.piie.a0.onestep.imp.var),3)

set6.result <- list(TMLE = data.frame(est = set6.TMLE.mean, var= set6.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set6.piie.a1.TMLE.mean, var= set6.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set6.piie.a0.TMLE.mean, var= set6.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set6.onestep.mean, var=set6.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set6.piie.a1.onestep.mean, var= set6.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set6.piie.a0.onestep.mean, var= set6.piie.a0.onestep.imp.var)
)


## Setting 7: M.cutoff2, Y.cutoff4 ====
set7.TMLE.est <- rep(NA,10); set7.TMLE.var <- rep(NA,10)
set7.onestep.est <- rep(NA,10); set7.onestep.var <- rep(NA,10)
set7.piie.a1.TMLE.est <- rep(NA,10); set7.piie.a1.TMLE.var <- rep(NA,10)
set7.piie.a0.TMLE.est <- rep(NA,10); set7.piie.a0.TMLE.var <- rep(NA,10)
set7.piie.a1.onestep.est <- rep(NA,10); set7.piie.a1.onestep.var <- rep(NA,10)
set7.piie.a0.onestep.est <- rep(NA,10); set7.piie.a0.onestep.var <- rep(NA,10)

set.seed(7)
for (m in 1:10){
  
  cat("Setting 7, iteration ", m, "\n")
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M.cutoff2"),
              outcome="Y.cutoff4", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
              lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
  
  set7.TMLE.est[m] <- est$TMLE$ATE
  set7.onestep.est[m] <- est$Onestep$ATE
  set7.piie.a1.TMLE.est[m] <- mean(dt.analysis$Y.cutoff4) - est$TMLE.Y1$estimated_psi
  set7.piie.a0.TMLE.est[m] <- mean(dt.analysis$Y.cutoff4) - est$TMLE.Y0$estimated_psi
  set7.piie.a1.onestep.est[m] <- mean(dt.analysis$Y.cutoff4) - est$Onestep.Y1$estimated_psi
  set7.piie.a0.onestep.est[m] <- mean(dt.analysis$Y.cutoff4) - est$Onestep.Y0$estimated_psi
  
  set7.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
  set7.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
  set7.piie.a1.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$TMLE.Y1$EIF}^2)/n
  set7.piie.a0.TMLE.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$TMLE.Y0$EIF}^2)/n
  set7.piie.a1.onestep.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$Onestep.Y1$EIF}^2)/n
  set7.piie.a0.onestep.var[m] <- mean({(dt.analysis$Y.cutoff4 - mean(dt.analysis$Y.cutoff4)) - est$Onestep.Y0$EIF}^2)/n
}

set7.TMLE.mean <- mean(set7.TMLE.est)
set7.onestep.mean <- mean(set7.onestep.est)
set7.piie.a1.TMLE.mean <- mean(set7.piie.a1.TMLE.est)
set7.piie.a0.TMLE.mean <- mean(set7.piie.a0.TMLE.est)
set7.piie.a1.onestep.mean <- mean(set7.piie.a1.onestep.est)
set7.piie.a0.onestep.mean <- mean(set7.piie.a0.onestep.est)

set7.TMLE.imp.var <- mean(set7.TMLE.var) + sum((set7.TMLE.est - set7.TMLE.mean)^2)/9
set7.onestep.imp.var <- mean(set7.onestep.var) + sum((set7.onestep.est - set7.onestep.mean)^2)/9
set7.piie.a1.TMLE.imp.var <- mean(set7.piie.a1.TMLE.var) + sum((set7.piie.a1.TMLE.est - set7.piie.a1.TMLE.mean)^2)/9
set7.piie.a0.TMLE.imp.var <- mean(set7.piie.a0.TMLE.var) + sum((set7.piie.a0.TMLE.est - set7.piie.a0.TMLE.mean)^2)/9
set7.piie.a1.onestep.imp.var <- mean(set7.piie.a1.onestep.var) + sum((set7.piie.a1.onestep.est - set7.piie.a1.onestep.mean)^2)/9
set7.piie.a0.onestep.imp.var <- mean(set7.piie.a0.onestep.var) + sum((set7.piie.a0.onestep.est - set7.piie.a0.onestep.mean)^2)/9

round(set7.TMLE.mean - 1.96*sqrt(set7.TMLE.imp.var),3); round(set7.TMLE.mean + 1.96*sqrt(set7.TMLE.imp.var),3)
round(set7.onestep.mean - 1.96*sqrt(set7.onestep.imp.var),3); round(set7.onestep.mean + 1.96*sqrt(set7.onestep.imp.var),3)
round(set7.piie.a1.TMLE.mean - 1.96*sqrt(set7.piie.a1.TMLE.imp.var),3); round(set7.piie.a1.TMLE.mean + 1.96*sqrt(set7.piie.a1.TMLE.imp.var),3)
round(set7.piie.a0.TMLE.mean - 1.96*sqrt(set7.piie.a0.TMLE.imp.var),3); round(set7.piie.a0.TMLE.mean + 1.96*sqrt(set7.piie.a0.TMLE.imp.var),3)
round(set7.piie.a1.onestep.mean - 1.96*sqrt(set7.piie.a1.onestep.imp.var),3); round(set7.piie.a1.onestep.mean + 1.96*sqrt(set7.piie.a1.onestep.imp.var),3)
round(set7.piie.a0.onestep.mean - 1.96*sqrt(set7.piie.a0.onestep.imp.var),3); round(set7.piie.a0.onestep.mean + 1.96*sqrt(set7.piie.a0.onestep.imp.var),3)

set7.result <- list(TMLE = data.frame(est = set7.TMLE.mean, var= set7.TMLE.imp.var), 
                    TMLE.piie.a1 = data.frame(est = set7.piie.a1.TMLE.mean, var= set7.piie.a1.TMLE.imp.var),
                    TMLE.piie.a0 = data.frame(est = set7.piie.a0.TMLE.mean, var= set7.piie.a0.TMLE.imp.var),
                    onestep = data.frame(est=set7.onestep.mean, var=set7.onestep.imp.var),
                    onestep.piie.a1 = data.frame(est = set7.piie.a1.onestep.mean, var= set7.piie.a1.onestep.imp.var),
                    onestep.piie.a0 = data.frame(est = set7.piie.a0.onestep.mean, var= set7.piie.a0.onestep.imp.var)
)


save(list = paste0('set', 1:7, '.result'), file="B_PROUD/estimation_v2.Rdata")