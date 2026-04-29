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
              lib = c("SL.bayesglm", "SL.earth" ,"SL.ranger", "SL.mean"), cvg.criteria = n^{-0.5})
  
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
       lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.mean"), cvg.criteria = n^{-0.5})
  
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

## organize the results====
library(huxtable)
library(dplyr)
options(scipen = 999)

load("B_PROUD/estimation_v2.Rdata")
load("output/estimation_v2.Rdata")

dat <- data.frame(Scenario= rep(c("Continuous",rep('Binary',6)),3),
                  M.cutoff=rep(c('-',rep('48 mins',3), rep('75 mins',3)),3),
                  Y.cutoff=rep(c('-',rep(2:4,2)),3),
                  Onestep=rep(NA,21),
                  TMLE=rep(NA,21)
                  )

for (i in 1:7){
  onestep.est <- get(paste0('set',i,'.result'))$onestep$est
  onestep.var <- get(paste0('set',i,'.result'))$onestep$var
  TMLE.est <- get(paste0('set',i,'.result'))$TMLE$est
  TMLE.var <- get(paste0('set',i,'.result'))$TMLE$var
  
  piie.a1.onestep.est <- get(paste0('set',i,'.result'))$onestep.piie.a1$est
  piie.a1.onestep.var <- get(paste0('set',i,'.result'))$onestep.piie.a1$var
  piie.a1.TMLE.est <- get(paste0('set',i,'.result'))$TMLE.piie.a1$est
  piie.a1.TMLE.var <- get(paste0('set',i,'.result'))$TMLE.piie.a1$var
  piie.a0.onestep.est <- get(paste0('set',i,'.result'))$onestep.piie.a0$est
  piie.a0.onestep.var <- get(paste0('set',i,'.result'))$onestep.piie.a0$var
  piie.a0.TMLE.est <- get(paste0('set',i,'.result'))$TMLE.piie.a0$est
  piie.a0.TMLE.var <- get(paste0('set',i,'.result'))$TMLE.piie.a0$var
  
  dat[i,4] <- paste0(round(onestep.est,3),", 95%CI (", round(onestep.est - 1.96*sqrt(onestep.var),3),", ", round(onestep.est + 1.96*sqrt(onestep.var),3),")")
  dat[i,5] <- paste0(round(TMLE.est,3),", 95%CI (", round(TMLE.est - 1.96*sqrt(TMLE.var),3),", ", round(TMLE.est + 1.96*sqrt(TMLE.var),3),")")
  
  dat[i+7,4] <- paste0(round(piie.a1.onestep.est,3),", 95%CI (", round(piie.a1.onestep.est - 1.96*sqrt(piie.a1.onestep.var),3),", ", round(piie.a1.onestep.est + 1.96*sqrt(piie.a1.onestep.var),3),")")
  dat[i+7,5] <- paste0(round(piie.a1.TMLE.est,3),", 95%CI (", round(piie.a1.TMLE.est - 1.96*sqrt(piie.a1.TMLE.var),3),", ", round(piie.a1.TMLE.est + 1.96*sqrt(piie.a1.TMLE.var),3),")")
  
  dat[i+14,4] <- paste0(round(piie.a0.onestep.est,3),", 95%CI (", round(piie.a0.onestep.est - 1.96*sqrt(piie.a0.onestep.var),3),", ", round(piie.a0.onestep.est + 1.96*sqrt(piie.a0.onestep.var),3),")")
  dat[i+14,5] <- paste0(round(piie.a0.TMLE.est,3),", 95%CI (", round(piie.a0.TMLE.est - 1.96*sqrt(piie.a0.TMLE.var),3),", ", round(piie.a0.TMLE.est + 1.96*sqrt(piie.a0.TMLE.var),3),")")
  
  
}

table1 <- as_hux(dat) %>%
  slice(-1) %>%  # Remove the first row
  insert_row("$M$ and $Y$ type","$M$ cutoff","$Y$ cutoff","One-step estimator","TMLE", after = 0) %>% # Insert first row
  insert_row("ATE","","","","", after = 1) %>% # Insert ATE row
  insert_row("PIIE $E(Y)-E(Y^{A,M^1})$","","","","", after = 9) %>% # Insert PIIE a1 row
  insert_row("PIIE $E(Y)-E(Y^{A,M^0})$","","","","", after = 17) %>% # Insert PIIE a1 row
  merge_cells(2, 1:5) %>% # Merge the first row
  merge_cells(10, 1:5) %>% # Merge the second row
  merge_cells(18, 1:5) %>% # Merge the third row
  set_escape_contents(1, 1:3, FALSE) %>%
  set_escape_contents(10, 1, FALSE) %>%
  set_escape_contents(18, 1, FALSE) %>%
  set_align(col=c(1,2,3), everywhere, "left") %>%
  set_align(col=c(4,5),everywhere,"center") %>%
  # set_align(everywhere,row=c(2,10,18),"center") %>%
  set_bold(everywhere,row=c(2,10,18)) %>% 
  set_valign(2:nrow(.), 1:5, "middle" ) %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere) %>%
  # set_italic(c(2,4,13,18),everywhere) %>%
  set_bottom_border(row = 1, col =1:ncol(.)) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "solid")) %>% 
  set_caption("One-step and TMLE estimates of the average causal effect of additional mobile stroke unit (MSU) care on modified Rankin scale (mRS)") %>%
  set_all_padding(1) %>% set_font_size(6)  %>% set_width(value=1) %>% set_col_width(c(0.15,0.1,0.1,0.325,0.325)) %>% set_wrap(T) %>% 
  insert_row(c("* Adopt one-step estimator $\\psi_1(\\hat{Q})$ and TMLE $\\psi_1(\\hat{Q}^\\star)$ under binary $M$, and adopt one-step estimator $\\psi_{2b}(\\hat{Q})$ and TMLE $\\psi_{2b}(\\hat{Q}^\\star)$ under continuous $M$.","","","",""), after = nrow(.)) %>% 
  merge_cells(nrow(.), 1:ncol(.)) %>% 
  set_escape_contents(nrow(.), 1:ncol(.), FALSE) %>% 
  set_bottom_border(row = nrow(.)-1,col = everywhere,brdr(1, "solid"))

table1

quick_latex(table1,file = "B_PROUD/table.tex")


#-----------------
dat <- data.frame(Scenario= c("Continuous",rep('Binary',6)),
                  M.cutoff=c('-',rep('48 mins',3), rep('75 mins',3)),
                  Y.cutoff=c('-',rep(2:4,2)),
                  Onestep=rep(NA,7),
                  TMLE=rep(NA,7)
)

for (i in 1:7){
  onestep.est <- get(paste0('set',i,'.result'))$onestep$est
  onestep.var <- get(paste0('set',i,'.result'))$onestep$var
  TMLE.est <- get(paste0('set',i,'.result'))$TMLE$est
  TMLE.var <- get(paste0('set',i,'.result'))$TMLE$var
  
  piie.a1.onestep.est <- get(paste0('set',i,'.result'))$onestep.piie.a1$est
  piie.a1.onestep.var <- get(paste0('set',i,'.result'))$onestep.piie.a1$var
  piie.a1.TMLE.est <- get(paste0('set',i,'.result'))$TMLE.piie.a1$est
  piie.a1.TMLE.var <- get(paste0('set',i,'.result'))$TMLE.piie.a1$var
  piie.a0.onestep.est <- get(paste0('set',i,'.result'))$onestep.piie.a0$est
  piie.a0.onestep.var <- get(paste0('set',i,'.result'))$onestep.piie.a0$var
  piie.a0.TMLE.est <- get(paste0('set',i,'.result'))$TMLE.piie.a0$est
  piie.a0.TMLE.var <- get(paste0('set',i,'.result'))$TMLE.piie.a0$var
  
  dat[i,4] <- paste0(round(onestep.est,3),", (", round(onestep.est - 1.96*sqrt(onestep.var),3),", ", round(onestep.est + 1.96*sqrt(onestep.var),3),")")
  dat[i,5] <- paste0(round(TMLE.est,3),", (", round(TMLE.est - 1.96*sqrt(TMLE.var),3),", ", round(TMLE.est + 1.96*sqrt(TMLE.var),3),")")
  
  dat[i,6] <- paste0(round(piie.a1.onestep.est,3),", (", round(piie.a1.onestep.est - 1.96*sqrt(piie.a1.onestep.var),3),", ", round(piie.a1.onestep.est + 1.96*sqrt(piie.a1.onestep.var),3),")")
  dat[i,7] <- paste0(round(piie.a1.TMLE.est,3),", (", round(piie.a1.TMLE.est - 1.96*sqrt(piie.a1.TMLE.var),3),", ", round(piie.a1.TMLE.est + 1.96*sqrt(piie.a1.TMLE.var),3),")")
  
  dat[i,8] <- paste0(round(piie.a0.onestep.est,3),", (", round(piie.a0.onestep.est - 1.96*sqrt(piie.a0.onestep.var),3),", ", round(piie.a0.onestep.est + 1.96*sqrt(piie.a0.onestep.var),3),")")
  dat[i,9] <- paste0(round(piie.a0.TMLE.est,3),", (", round(piie.a0.TMLE.est - 1.96*sqrt(piie.a0.TMLE.var),3),", ", round(piie.a0.TMLE.est + 1.96*sqrt(piie.a0.TMLE.var),3),")")
  
  
}

table1 <- as_hux(dat) %>%
  slice(-1) %>%  # Remove the first row
  insert_row("$M$ and $Y$ type","Cutoffs","","ATE","","PIIE at $a_0=1$","","PIIE at $a_0=0$","", after = 0) %>% # Insert first row
  insert_row("","$M$","$Y$","One-step","TMLE","One-step","TMLE","One-step","TMLE", after = 1) %>% # Insert first row
  set_escape_contents(1:2, everywhere, FALSE) %>%
  merge_cells(1, 2:3) %>% 
  merge_cells(1, 4:5) %>% 
  merge_cells(1, 6:7) %>% 
  merge_cells(1, 8:9) %>%
  set_align(col=c(1), everywhere, "left") %>%
  set_align(col=c(2,3,4,5,6,7,8,9),everywhere,"center") %>%
  set_valign(2:nrow(.), 1:9, "middle" ) %>%
  set_tb_padding(1, everywhere, 0) %>% 
  set_bold(1, everywhere) %>%
  # set_italic(c(2,4,13,18),everywhere) %>%
  set_bottom_border(row = 1, col =1:ncol(.)) %>% 
  set_bottom_border(row=2,col=2:ncol(.)) %>% 
  set_top_border(row=1,col=everywhere,brdr(1, "solid")) %>% 
  set_caption("One-step and TMLE estimates of the average causal effect of additional mobile stroke unit (MSU) care on modified Rankin scale (mRS)") %>%
  set_all_padding(1) %>% set_font_size(6)  %>% set_width(value=1) %>% 
  insert_row(c("* Adopt one-step estimator $\\psi_1(\\hat{Q})$ and TMLE $\\psi_1(\\hat{Q}^\\star)$ under binary $M$, and adopt one-step estimator $\\psi_{2b}(\\hat{Q})$ and TMLE $\\psi_{2b}(\\hat{Q}^\\star)$ under continuous $M$. Estimates are reported, together with $95\\%$ confidence interval.","","","","","","","",""), after = nrow(.)) %>% 
  merge_cells(nrow(.), 1:ncol(.)) %>% 
  set_escape_contents(nrow(.), 1:ncol(.), FALSE) %>% 
  set_bottom_border(row = nrow(.)-1,col = everywhere,brdr(1, "solid"))

table1

quick_latex(table1,file = "B_PROUD/table.tex")


# # Use original continuous M and treat Y as continuous================
# 
# ## compare whether the original data (dt2) with continuous M agree with the organized categorized data (dt1)
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
# 
# set1.TMLE.est <- rep(NA,10); set1.TMLE.var <- rep(NA,10)
# set1.onestep.est <- rep(NA,10); set1.onestep.var <- rep(NA,10)
# 
# set.seed(7)
# for (m in 1:10){
#   
#   dt.analysis <- dt %>% filter(.imp==m)
#   n <- nrow(dt.analysis)
#   
#   est <- estfd(a=c(1,0),data=dt.analysis,treatment="A", mediators=c("M1","M2"),
#               outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
#               lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = n^{-0.5})
#   
#   set1.TMLE.est[m] <- est$TMLE$ATE
#   set1.onestep.est[m] <- est$Onestep$ATE
#   
#   set1.TMLE.var[m] <- mean(est$TMLE$EIF^2)/n
#   set1.onestep.var[m] <- mean(est$Onestep$EIF^2)/n
# }
# 
# set1.TMLE.mean <- mean(set1.TMLE.est)
# set1.onestep.mean <- mean(set1.onestep.est)
# set1.TMLE.imp.var <- mean(set1.TMLE.var) + sum((set1.TMLE.est - set1.TMLE.mean)^2)/9
# set1.onestep.imp.var <- mean(set1.onestep.var) + sum((set1.onestep.est - set1.onestep.mean)^2)/9
# 
# round(set1.TMLE.mean - 1.96*sqrt(set1.TMLE.imp.var),3); round(set1.TMLE.mean + 1.96*sqrt(set1.TMLE.imp.var),3)
# round(set1.onestep.mean - 1.96*sqrt(set1.onestep.imp.var),3); round(set1.onestep.mean + 1.96*sqrt(set1.onestep.imp.var),3)
# 
# set_contiM_contiY.result <- list(TMLE = data.frame(est = set1.TMLE.mean, var= set1.TMLE.imp.var), onestep = data.frame(est=set1.onestep.mean, var=set1.onestep.imp.var))
# save(list = c('set_contiM_contiY.result'), file="B_PROUD/estimation_contiM_contiY.Rdata")


# Use original continuous M and ordinal Y================

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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  #a=0
  est_Y0_a0 <- estfd(a=0,data=dt_Y0,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger","SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y0)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y1_a0 <- estfd(a=0,data=dt_Y1,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y1)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y2_a0 <- estfd(a=0,data=dt_Y2,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y2)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y3_a0 <- estfd(a=0,data=dt_Y3,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y3)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y4_a0 <- estfd(a=0,data=dt_Y4,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y4)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
  #a=0
  est_Y5_a0 <- estfd(a=0,data=dt_Y5,treatment="A", mediators=c("M1","M2"),
                    outcome="Y", covariates=c("X1","X2"), estimator=c('onestep','tmle'), mediator.method="bayes",superlearner = T, crossfit = T,K=5,
                    lib = c("SL.glm", "SL.earth" ,"SL.ranger", "SL.xgboost", "SL.mean"), cvg.criteria = nrow(dt_Y5)^{-0.5}, truncate_lower = 0.01, truncate_upper = 0.99)
  
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


set1.TMLE.dim.est <- rep(NA,10); set1.TMLE.dim.var <- rep(NA,10)
set1.TMLE.piie.a1.est <- rep(NA,10); set1.TMLE.piie.a1.var <- rep(NA,10)
set1.TMLE.piie.a0.est <- rep(NA,10); set1.TMLE.piie.a0.var <- rep(NA,10)
set1.TMLE.lor.est <- rep(NA,10); set1.TMLE.lor.var <- rep(NA,10)
set1.onestep.dim.est <- rep(NA,10); set1.onestep.dim.var <- rep(NA,10)
set1.onestep.piie.a1.est <- rep(NA,10); set1.onestep.piie.a1.var <- rep(NA,10)
set1.onestep.piie.a0.est <- rep(NA,10); set1.onestep.piie.a0.var <- rep(NA,10)
set1.onestep.lor.est <- rep(NA,10); set1.onestep.lor.var <- rep(NA,10)

TMLE.PMF.Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
TMLE.PMF.Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")

set.seed(7)
for (m in 1:10){
  
  cat('Imputation', m, 'of 10\n')
  
  dt.analysis <- dt %>% filter(.imp==m)
  n <- nrow(dt.analysis)
  
  # bootstrap for DIM and LOR
  estimation.m <- TMLE_m(dt.analysis)
  
  set1.TMLE.dim.est[m] <- estimation.m$dim_tmle
  set1.TMLE.piie.a1.est[m] <- estimation.m$piie_tmle_a1
  set1.TMLE.piie.a0.est[m] <- estimation.m$piie_tmle_a0
  set1.TMLE.lor.est[m] <- estimation.m$lor_tmle
  
  set1.onestep.dim.est[m] <- estimation.m$dim_onestep
  set1.onestep.piie.a1.est[m] <- estimation.m$piie_onestep_a1
  set1.onestep.piie.a0.est[m] <- estimation.m$piie_onestep_a0
  set1.onestep.lor.est[m] <- estimation.m$lor_onestep
  
  set1.TMLE.dim.var[m] <- mean(estimation.m$dim_tmle_EIF^2)/nrow(dt.analysis)
  set1.TMLE.piie.a1.var[m] <- mean(estimation.m$piie_tmle_a1_EIF^2)/nrow(dt.analysis)
  set1.TMLE.piie.a0.var[m] <- mean(estimation.m$piie_tmle_a0_EIF^2)/nrow(dt.analysis)
  
  set1.onestep.dim.var[m] <-mean(estimation.m$dim_onestep_EIF^2)/nrow(dt.analysis)
  set1.onestep.piie.a1.var[m] <- mean(estimation.m$piie_onestep_a1_EIF^2)/nrow(dt.analysis)
  set1.onestep.piie.a0.var[m] <- mean(estimation.m$piie_onestep_a0_EIF^2)/nrow(dt.analysis)
  
  TMLE.PMF.Ya1[m,] <- c(m, estimation.m$tmle_Y0_a1 ,estimation.m$tmle_Y1_a1, estimation.m$tmle_Y2_a1, estimation.m$tmle_Y3_a1, estimation.m$tmle_Y4_a1, estimation.m$tmle_Y5_a1, estimation.m$tmle_Y6_a1)
  TMLE.PMF.Ya0[m,] <- c(m, estimation.m$tmle_Y0_a0 ,estimation.m$tmle_Y1_a0, estimation.m$tmle_Y2_a0, estimation.m$tmle_Y3_a0, estimation.m$tmle_Y4_a0, estimation.m$tmle_Y5_a0, estimation.m$tmle_Y6_a0)
  onestep.PMF.Ya1[m,] <- c(m, estimation.m$onestep_Y0_a1 ,estimation.m$onestep_Y1_a1, estimation.m$onestep_Y2_a1, estimation.m$onestep_Y3_a1, estimation.m$onestep_Y4_a1, estimation.m$onestep_Y5_a1, estimation.m$onestep_Y6_a1)
  onestep.PMF.Ya0[m,] <- c(m, estimation.m$onestep_Y0_a0 ,estimation.m$onestep_Y1_a0, estimation.m$onestep_Y2_a0, estimation.m$onestep_Y3_a0, estimation.m$onestep_Y4_a0, estimation.m$onestep_Y5_a0, estimation.m$onestep_Y6_a0)
  
}

set1.TMLE.dim.mean <- mean(set1.TMLE.dim.est)
set1.onestep.dim.mean <- mean(set1.onestep.dim.est)

set1.TMLE.piie.a1.mean <- mean(set1.TMLE.piie.a1.est)
set1.TMLE.piie.a0.mean <- mean(set1.TMLE.piie.a0.est)
set1.onestep.piie.a1.mean <- mean(set1.onestep.piie.a1.est)
set1.onestep.piie.a0.mean <- mean(set1.onestep.piie.a0.est)

set1.TMLE.lor.mean <- mean(set1.TMLE.lor.est)
set1.onestep.lor.mean <- mean(set1.onestep.lor.est)

set1.TMLE.dim.imp.var <- mean(set1.TMLE.dim.var) + sum((set1.TMLE.dim.est - set1.TMLE.dim.mean)^2)/9
set1.onestep.dim.imp.var <- mean(set1.onestep.dim.var) + sum((set1.onestep.dim.est - set1.onestep.dim.mean)^2)/9
set1.TMLE.piie.a1.imp.var <- mean(set1.TMLE.piie.a1.var) + sum((set1.TMLE.piie.a1.est - set1.TMLE.piie.a1.mean)^2)/9
set1.TMLE.piie.a0.imp.var <- mean(set1.TMLE.piie.a0.var) + sum((set1.TMLE.piie.a0.est - set1.TMLE.piie.a0.mean)^2)/9
set1.onestep.piie.a1.imp.var <- mean(set1.onestep.piie.a1.var) + sum((set1.onestep.piie.a1.est - set1.onestep.piie.a1.mean)^2)/9
set1.onestep.piie.a0.imp.var <- mean(set1.onestep.piie.a0.var) + sum((set1.onestep.piie.a0.est - set1.onestep.piie.a0.mean)^2)/9


round(set1.TMLE.dim.mean - 1.96*sqrt(set1.TMLE.dim.imp.var),3); round(set1.TMLE.dim.mean + 1.96*sqrt(set1.TMLE.dim.imp.var),3)
round(set1.onestep.dim.mean - 1.96*sqrt(set1.onestep.dim.imp.var),3); round(set1.onestep.dim.mean + 1.96*sqrt(set1.onestep.dim.imp.var),3)


set_contiM_ordinalY.result <- list(TMLE = data.frame(dim = set1.TMLE.dim.mean, dim.var= set1.TMLE.dim.imp.var, lor = set1.TMLE.lor.mean, 
                                                     piie.a1= set1.TMLE.piie.a1.mean, piie.a0 = set1.TMLE.piie.a0.mean,
                                                     piie.a1.var = set1.TMLE.piie.a1.imp.var, piie.a0.var = set1.TMLE.piie.a0.imp.var),
                                   onestep = data.frame(dim = set1.onestep.dim.mean, dim.var= set1.onestep.dim.imp.var, lor = set1.onestep.lor.mean,
                                                        piie.a1= set1.onestep.piie.a1.mean, piie.a0 = set1.onestep.piie.a0.mean,
                                                     piie.a1.var = set1.onestep.piie.a1.imp.var, piie.a0.var = set1.onestep.piie.a0.imp.var))
save(list = c('set_contiM_ordinalY.result',
              'set1.TMLE.dim.mean','set1.onestep.dim.mean',
              'set1.TMLE.piie.a1.mean','set1.TMLE.piie.a0.mean',
              'set1.TMLE.piie.a1.imp.var','set1.TMLE.piie.a0.imp.var',
              'set1.onestep.piie.a1.mean','set1.onestep.piie.a0.mean',
              'set1.onestep.piie.a1.imp.var','set1.onestep.piie.a0.imp.var',
              'set1.TMLE.lor.mean','set1.onestep.lor.mean',
              'set1.TMLE.dim.imp.var','set1.onestep.dim.imp.var',
              'TMLE.PMF.Ya1','TMLE.PMF.Ya0',
              'onestep.PMF.Ya1','onestep.PMF.Ya0'), file="B_PROUD/estimation_contiM_ordinalY.Rdata")
