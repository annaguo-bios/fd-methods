EYa.binaryZ.binaryM <- function(a,data,treatment, mediator, outcome, anchor,
                        onestep=TRUE, superlearner=TRUE,crossfit=FALSE,K=5,
                        lib = c("SL.glm","SL.earth","SL.ranger","SL.mean"), n.iter=500, cvg.criteria=0.01,
                        formulaY="Y ~ .", formulaA="A ~ .", formulaM="M~.", linkY_binary="logit", linkA="logit", linkM_binary="logit",
                        truncate_lower=0.01, truncate_upper=0.99,
                        verbose=T){

  # attach(data, warn.conflicts=FALSE)

  n <- nrow(data)
  # Variables
  A <- data[,treatment]
  M <- data[,mediator]
  Z <- data[,anchor]
  Y <- data[,outcome]
  
  # p(Z)
  density.z <-function(z){
    if(is.null(z)){NA}else{mean(Z==z)}
  }
  
  # p(A|Z)
  density.a.z <- function(a,z){
    
    p.a1.z1 <- mean(A[Z==1]==1)
    p.a1.z0 <- mean(A[Z==0]==1)
    
    out <- dbinom(a, 1, (z*p.a1.z1 + (1-z)*p.a1.z0))
    
    return(out)
  }
  
  # p(M|A,Z)
  
  density.m.az <- function(m,a,z){
    
    p.m1.a1z1 <- mean(M[A==1 & Z==1]==1)
    p.m1.a1z0 <- mean(M[A==1 & Z==0]==1)
    p.m1.a0z1 <- mean(M[A==0 & Z==1]==1)
    p.m1.a0z0 <- mean(M[A==0 & Z==0]==1)
    
    out <- dbinom(m, 1, {a*z*p.m1.a1z1 + a*(1-z)*p.m1.a1z0 + (1-a)*z*p.m1.a0z1 + (1-z)*(1-a)*p.m1.a0z0})
    
    out[out==0] <- 0.001
    
    return(out)
    
  }
  
  # E(Y|M,A,Z)
  safe_mean <- function(x) {
    if (length(x) == 0) return(0)
    mean(x)
  }
  E.y.maz <- function(m,a,z){
    
    Ey.m1a1z1 <- safe_mean(Y[M==1 & A==1 & Z==1])
    Ey.m1a1z0 <- safe_mean(Y[M==1 & A==1 & Z==0])
    Ey.m1a0z1 <- safe_mean(Y[M==1 & A==0 & Z==1])
    Ey.m1a0z0 <- safe_mean(Y[M==1 & A==0 & Z==0])
    
    Ey.m0a1z1 <- safe_mean(Y[M==0 & A==1 & Z==1])
    Ey.m0a1z0 <- safe_mean(Y[M==0 & A==1 & Z==0])
    Ey.m0a0z1 <- safe_mean(Y[M==0 & A==0 & Z==1])
    Ey.m0a0z0 <- safe_mean(Y[M==0 & A==0 & Z==0])
    
    out <- m*{a*z*Ey.m1a1z1 + a*(1-z)*Ey.m1a1z0 + (1-a)*z*Ey.m1a0z1 + (1-a)*(1-z)*Ey.m1a0z0} + (1-m)*{a*z*Ey.m0a1z1 + a*(1-z)*Ey.m0a1z0 + (1-a)*z*Ey.m0a0z1 + (1-a)*(1-z)*Ey.m0a0z0}
    
    return(out)
    
  }
  
  ##################################################################
  #################### One-step estimator ##########################
  ##################################################################
  
  f.onestep.binaryZM <- function(a,z=NULL){ # 1. if z=null, then return estimand from the regular ID functional where z is not fixed, 2. if z is not null, then return the estimand from the ID functional where z is fixed
    
    # outcome regression E[Y|M,A,z]
    if (!is.null(z)){
    mu.a1m1_z <- E.y.maz(1,1,z) # A=1, M=1, Z=z
    mu.a1m0_z <- E.y.maz(0,1,z) # A=1, M=0, Z=z
    mu.a0m1_z <- E.y.maz(1,0,z) # A=0, M=1, Z=z
    mu.a0m0_z <- E.y.maz(0,0,z) # A=0, M=0, Z=z
    mu_z <- E.y.maz(M,A,z) # A=A, M=M, Z=z
    }
    
    mu <- E.y.maz(M,A,Z) # A=A, M=M, Z=Z
    mu.a1m1 <- E.y.maz(1,1,Z) # A=1, M=1, Z=Z
    mu.a1m0 <- E.y.maz(0,1,Z) # A=1, M=0, Z=Z
    mu.a0m1 <- E.y.maz(1,0,Z) # A=0, M=1, Z=Z
    mu.a0m0 <- E.y.maz(0,0,Z) # A=0, M=0, Z=Z
    
    # mediator density p(M|A=a,Z)
    p.m1.az1 <- density.m.az(m=1,a,z=1) #p(M=1|a,Z=1)
    p.m1.az0 <- density.m.az(m=1,a,z=0) #p(M=1|a,Z=0)
    p.m0.az1 <- density.m.az(m=0,a,z=1) #p(M=0|a,Z=1)
    p.m0.az0 <- density.m.az(m=0,a,z=0) #p(M=0|a,Z=0
    
    p.M.aZ <- density.m.az(M,a,Z) # p(M|A=a,Z)
    p.M.AZ <- density.m.az(M,A,Z) # p(M|A,Z)
    p.m1.aZ <- density.m.az(m=1,a,Z) # p(M=1|A=a,Z)
    p.m0.aZ <- density.m.az(m=0,a,Z) # p(M=0|A=a,Z)
    p.M.az1 <- density.m.az(M,a,z=1) # p(M|A=a,Z=1)
    p.M.az0 <- density.m.az(M,a,z=0) # p(M|A=a,Z=0)
    if (!is.null(z)){
      p.M.az <- density.m.az(M,a,z) # p(M|A=a,Z=z)
      p.M.Az <- density.m.az(M,A,z) # p(M|A,Z=z)
      } 
    
    # propensity score p(A|Z=z)
    if (!is.null(z)){
    p.a1.z <- density.a.z(a=1,z) # p(a=1|Z=z)
    p.a0.z <- density.a.z(a=0,z) # p(a=0|Z=z)
    }
    
    p.a.Z <- density.a.z(a,Z) # p(A=a|Z)
    p.a1.Z <- density.a.z(a=1,Z) # p(A=1|Z)
    p.a0.Z <- density.a.z(a=0,Z) # p(A=1|Z)
    
    # Z density
    p.z1 <- density.z(1) # p(Z=1)
    p.z0 <- density.z(0) # p(Z=0)
    if (!is.null(z)){p.z <- density.z(z)} # p(Z=z)
    
    # int E(Y|M,A,z)p(M|A=a,Z)p(A|z)p(Z)dMdAdZ
    if (is.null(z)){
      
      theta_z <- (mu.a1m1)*p.m1.aZ*p.a1.Z+ # M=1,A=1
        (mu.a0m1)*p.m1.aZ*p.a0.Z+ # M=1, A=0
        (mu.a1m0)*p.m0.aZ*p.a1.Z+ # M=0, A=1
        (mu.a0m0)*p.m0.aZ*p.a0.Z # M=0, A=0 
      
      psi <- mean(theta_z)
      
      # true EIF 
      EIF.Y <- {p.M.aZ/p.M.AZ}*(Y-mu)
      EIF.M <- I(A==a)/p.a.Z*(mu.a1m1*p.a1.Z+mu.a0m1*p.a0.Z-mu.a1m0*p.a1.Z-mu.a0m0*p.a0.Z)*(M-p.m1.aZ)
      EIF.A <- (mu.a1m1*p.m1.aZ+mu.a1m0*p.m0.aZ-mu.a0m1*p.m1.aZ-mu.a0m0*p.m0.aZ)*(A-p.a1.Z)
      EIF.Z <- theta_z - psi
      
      EIF <- EIF.Y+EIF.M+EIF.A+EIF.Z
      
      psi <- mean(EIF.Y+EIF.M+EIF.A) + mean(theta_z)
      
      
    }else{
      
      theta_z <- (mu.a1m1_z)*p.m1.aZ*p.a1.z+ # M=1,A=1
        (mu.a0m1_z)*p.m1.aZ*p.a0.z+ # M=1, A=0
        (mu.a1m0_z)*p.m0.aZ*p.a1.z+ # M=0, A=1
        (mu.a0m0_z)*p.m0.aZ*p.a0.z # M=0, A=0 
      
      psi <- mean(theta_z)
      
      # true EIF
      EIF.Y <- I(Z==z)*{p.M.az1*p.z1 + p.M.az0*p.z0}/{p.M.Az*p.z}*(Y-mu_z)
      EIF.M <- I(A==a)/p.a.Z*(mu.a1m1_z*p.a1.z+mu.a0m1_z*p.a0.z-mu.a1m0_z*p.a1.z-mu.a0m0_z*p.a0.z)*(M-p.m1.aZ)
      EIF.A <- I(Z==z)/p.z*{(mu.a1m1_z*p.m1.az1+mu.a1m0_z*p.m0.az1-mu.a0m1_z*p.m1.az1-mu.a0m0_z*p.m0.az1)*p.z1+(mu.a1m1_z*p.m1.az0+mu.a1m0_z*p.m0.az0-mu.a0m1_z*p.m1.az0-mu.a0m0_z*p.m0.az0)*p.z0}*(A-p.a1.z)
      EIF.Z <- theta_z - psi
      
      EIF <- EIF.Y + EIF.M + EIF.A + EIF.Z
      
      psi <- mean(EIF.Y + EIF.M + EIF.A) + mean(theta_z)
      
    }
    
    # confidence interval
    lower.ci <- psi-1.96*sqrt(mean(EIF^2)/n)
    upper.ci <- psi+1.96*sqrt(mean(EIF^2)/n)
    
    
    return(list(EIF=EIF, estimated=psi, lower.ci=lower.ci, upper.ci=upper.ci))
    
  } # end of f.EIF.binary function
  
  
  ## calculate the truth under different levels of z
  out_z1 <- f.onestep.binaryZM(a,1)
  out_z0 <- f.onestep.binaryZM(a,0)
  out_allz <- f.onestep.binaryZM(a)
  
  estimated_psi_z1 <- out_z1$estimated
  estimated_psi_z0 <- out_z0$estimated
  estimated_psi_allz <- out_allz$estimated
  
  # EIF under different Z
  EIF_z1 <- out_z1$EIF
  EIF_z0 <- out_z0$EIF
  EIF_allZ <- out_allz$EIF
  
  if(verbose){print("Z is binary. Computing one-step estimators at Z=1, Z=0, and Z=Z.")}
  
  onestep.out <- list(out.z1=out_z1, out.z0=out_z0, out.all.z=out_allz)
  
  
  ## calculate the truth under different levels of z
  out_z1 <- f.onestep.binaryZM(a,1)
  out_z0 <- f.onestep.binaryZM(a,0)
  out_allz <- f.onestep.binaryZM(a)
  
  estimated_psi_z1 <- out_z1$estimated
  estimated_psi_z0 <- out_z0$estimated
  estimated_psi_allz <- out_allz$estimated
  
  # EIF under different Z
  EIF_z1 <- out_z1$EIF
  EIF_z0 <- out_z0$EIF
  EIF_allZ <- out_allz$EIF
  
  if(verbose){print("Z is binary. Computing one-step estimators at Z=1, Z=0, and Z=Z.")}
  
  onestep.out <- list(out.z1=out_z1, out.z0=out_z0, out.all.z=out_allz)
  
  tmle.out <- list(out.z1=out_z1, out.z0=out_z0, out.all.z=out_allz)
  
  out <- list(Onestep=onestep.out, TMLE=tmle.out)
  
  return(out)
  
}
