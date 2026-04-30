
TMLE.dim.est <- rep(NA,10); TMLE.dim.var <- rep(NA,10)
TMLE.ATT.est <- rep(NA,10); TMLE.ATT.var <- rep(NA,10)
TMLE.piie.a1.est <- rep(NA,10); TMLE.piie.a1.var <- rep(NA,10)
TMLE.piie.a0.est <- rep(NA,10); TMLE.piie.a0.var <- rep(NA,10)

onestep.dim.est <- rep(NA,10); onestep.dim.var <- rep(NA,10)
onestep.ATT.est <- rep(NA,10); onestep.ATT.var <- rep(NA,10)
onestep.piie.a1.est <- rep(NA,10); onestep.piie.a1.var <- rep(NA,10)
onestep.piie.a0.est <- rep(NA,10); onestep.piie.a0.var <- rep(NA,10)

TMLE.PMF.Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
TMLE.PMF.Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
TMLE.PMF.ATT_Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.ATT_Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
TMLE.PMF.ATT_Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(TMLE.PMF.ATT_Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.ATT_Ya1 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.ATT_Ya1) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")
onestep.PMF.ATT_Ya0 <- data.frame(matrix(ncol = 8, nrow = 10)); colnames(onestep.PMF.ATT_Ya0) <- c("m", "E.Y0", "E.Y1", "E.Y2", "E.Y3", "E.Y4", "E.Y5", "E.Y6")

set.seed(7)
for (m in 1:10){
  
  load(paste0('output/output_', m, '.RData'))
  
  TMLE.dim.est[m] <- TMLE.dim.est_m
  TMLE.ATT.est[m] <- TMLE.ATT.est_m
  TMLE.piie.a1.est[m] <- TMLE.piie.a1.est_m
  TMLE.piie.a0.est[m] <- TMLE.piie.a0.est_m

  
  onestep.dim.est[m] <- onestep.dim.est_m
  onestep.ATT.est[m] <- onestep.ATT.est_m
  onestep.piie.a1.est[m] <- onestep.piie.a1.est_m
  onestep.piie.a0.est[m] <- onestep.piie.a0.est_m

  
  TMLE.dim.var[m] <- TMLE.dim.var_m
  TMLE.ATT.var[m] <- TMLE.ATT.var_m
  TMLE.piie.a1.var[m] <- TMLE.piie.a1.var_m
  TMLE.piie.a0.var[m] <- TMLE.piie.a0.var_m
  
  onestep.dim.var[m] <- onestep.dim.var_m
  onestep.ATT.var[m] <- onestep.ATT.var_m
  onestep.piie.a1.var[m] <- onestep.piie.a1.var_m
  onestep.piie.a0.var[m] <- onestep.piie.a0.var_m
  
  TMLE.PMF.Ya1[m,] <- TMLE.PMF.Ya1_m
  TMLE.PMF.Ya0[m,] <- TMLE.PMF.Ya0_m
  onestep.PMF.Ya1[m,] <- onestep.PMF.Ya1_m
  onestep.PMF.Ya0[m,] <- onestep.PMF.Ya0_m
  TMLE.PMF.ATT_Ya1[m,] <- TMLE.PMF.ATT_Ya1_m
  TMLE.PMF.ATT_Ya0[m,] <- TMLE.PMF.ATT_Ya0_m
  onestep.PMF.ATT_Ya1[m,] <- onestep.PMF.ATT_Ya1_m
  onestep.PMF.ATT_Ya0[m,] <- onestep.PMF.ATT_Ya0_m
  
}

TMLE.dim.mean <- mean(TMLE.dim.est)
onestep.dim.mean <- mean(onestep.dim.est)

TMLE.ATT.mean <- mean(TMLE.ATT.est)
onestep.ATT.mean <- mean(onestep.ATT.est)

TMLE.piie.a1.mean <- mean(TMLE.piie.a1.est)
TMLE.piie.a0.mean <- mean(TMLE.piie.a0.est)
onestep.piie.a1.mean <- mean(onestep.piie.a1.est)
onestep.piie.a0.mean <- mean(onestep.piie.a0.est)

TMLE.dim.imp.var <- mean(TMLE.dim.var) + sum((TMLE.dim.est - TMLE.dim.mean)^2)/9
onestep.dim.imp.var <- mean(onestep.dim.var) + sum((onestep.dim.est - onestep.dim.mean)^2)/9
TMLE.ATT.imp.var <- mean(TMLE.ATT.var) + sum((TMLE.ATT.est - TMLE.ATT.mean)^2)/9
onestep.ATT.imp.var <- mean(onestep.ATT.var) + sum((onestep.ATT.est - onestep.ATT.mean)^2)/9
TMLE.piie.a1.imp.var <- mean(TMLE.piie.a1.var) + sum((TMLE.piie.a1.est - TMLE.piie.a1.mean)^2)/9
TMLE.piie.a0.imp.var <- mean(TMLE.piie.a0.var) + sum((TMLE.piie.a0.est - TMLE.piie.a0.mean)^2)/9
onestep.piie.a1.imp.var <- mean(onestep.piie.a1.var) + sum((onestep.piie.a1.est - onestep.piie.a1.mean)^2)/9
onestep.piie.a0.imp.var <- mean(onestep.piie.a0.var) + sum((onestep.piie.a0.est - onestep.piie.a0.mean)^2)/9


set_contiM_ordinalY.result <- list(TMLE = data.frame(dim = TMLE.dim.mean, dim.var= TMLE.dim.imp.var,
                                                     ATT = TMLE.ATT.mean, ATT.var = TMLE.ATT.imp.var, 
                                                     piie.a1= TMLE.piie.a1.mean, piie.a0 = TMLE.piie.a0.mean,
                                                     piie.a1.var = TMLE.piie.a1.imp.var, piie.a0.var = TMLE.piie.a0.imp.var),
                                   onestep = data.frame(dim = onestep.dim.mean, dim.var= onestep.dim.imp.var,
                                                        ATT = onestep.ATT.mean, ATT.var = onestep.ATT.imp.var,
                                                        piie.a1= onestep.piie.a1.mean, piie.a0 = onestep.piie.a0.mean,
                                                        piie.a1.var = onestep.piie.a1.imp.var, piie.a0.var = onestep.piie.a0.imp.var))

round(colMeans(TMLE.PMF.Ya1),2)*100


# ATE #
# TMLE
cat('ATE TMLE:', round(TMLE.dim.mean,3),' (95\\% CI: (', round(TMLE.dim.mean-1.96*sqrt(TMLE.dim.imp.var),3),',', round(TMLE.dim.mean+1.96*sqrt(TMLE.dim.imp.var),3),'))\n')


# Onestep
cat('ATE onestep:',round(onestep.dim.mean,3),' (95\\% CI: (', round(onestep.dim.mean-1.96*sqrt(onestep.dim.imp.var),3),',', round(onestep.dim.mean+1.96*sqrt(onestep.dim.imp.var),3),'))\n',sep = '')

# ATT #
# TMLE
cat('ATT TMLE',round(TMLE.ATT.mean,3),' (95\\% CI: (', round(TMLE.ATT.mean-1.96*sqrt(TMLE.ATT.imp.var),3),',', round(TMLE.ATT.mean+1.96*sqrt(TMLE.ATT.imp.var),3),'))\n',sep = '')

# Onestep
cat('ATT onestep',round(onestep.ATT.mean,3),' (95\\% CI: (', round(onestep.ATT.mean-1.96*sqrt(onestep.ATT.imp.var),3),',', round(onestep.ATT.mean+1.96*sqrt(onestep.ATT.imp.var),3),'))\n',sep = '')

# piie a1 #
# TMLE
cat('PIIE a1 TMLE:',round(TMLE.piie.a1.mean,3),' (95\\% CI: (', round(TMLE.piie.a1.mean-1.96*sqrt(TMLE.piie.a1.imp.var),3),',', round(TMLE.piie.a1.mean+1.96*sqrt(TMLE.piie.a1.imp.var),3),')\n',sep = '')

# Onestep
cat('PIIE a1 onestep:',round(onestep.piie.a1.mean,3),' (95\\% CI: (', round(onestep.piie.a1.mean-1.96*sqrt(onestep.piie.a1.imp.var),3),',', round(onestep.piie.a1.mean+1.96*sqrt(onestep.piie.a1.imp.var),3),')\n',sep = '')

# piie a0 #
# TMLE
cat('PIIE a0 TMLE:',round(TMLE.piie.a0.mean,3),' (95\\% CI: (', round(TMLE.piie.a0.mean-1.96*sqrt(TMLE.piie.a0.imp.var),3),',', round(TMLE.piie.a0.mean+1.96*sqrt(TMLE.piie.a0.imp.var),3),')\n',sep = '')

# Onestep
cat('PIIE a0 onestep:',round(onestep.piie.a0.mean,3),' (95\\% CI: (', round(onestep.piie.a0.mean-1.96*sqrt(onestep.piie.a0.imp.var),3),',', round(onestep.piie.a0.mean+1.96*sqrt(onestep.piie.a0.imp.var),3),')\n',sep = '')


save(list = c('set_contiM_ordinalY.result',
              'TMLE.dim.mean','onestep.dim.mean',
              'TMLE.ATT.mean','onestep.ATT.mean',
              'TMLE.piie.a1.mean','TMLE.piie.a0.mean',
              'TMLE.piie.a1.imp.var','TMLE.piie.a0.imp.var',
              'onestep.piie.a1.mean','onestep.piie.a0.mean',
              'onestep.piie.a1.imp.var','onestep.piie.a0.imp.var',
              'TMLE.dim.imp.var','onestep.dim.imp.var',
              'TMLE.ATT.imp.var','onestep.ATT.imp.var',
              'TMLE.PMF.Ya1','TMLE.PMF.Ya0',
              'onestep.PMF.Ya1','onestep.PMF.Ya0',
              'TMLE.PMF.ATT_Ya1','TMLE.PMF.ATT_Ya0',
              'onestep.PMF.ATT_Ya1','onestep.PMF.ATT_Ya0'
              ), file="output/estimation_contiM_ordinalY.Rdata")
