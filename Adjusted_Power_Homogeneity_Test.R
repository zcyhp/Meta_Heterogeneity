#specify the packages of interest
packages = c("Metrics","lme4", "metafor","BiasedUrn", "MCMCpack", "aod", "CompQuadForm")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
###########################################################
# Function to generate n, number of subject in each group #
###########################################################
nct_gen = function(low, up, R, K, M){
  nc = matrix(NA, M, K)
  nt = matrix(NA, M, K)
  nc = matrix(runif(K*M, min=low, max=up), M, K)
  nt = matrix(2^(rnorm(K*M, log2(R), sqrt(0.5))), M, K)*nc
  return(list(round(nc),round(nt)))
}
#########################################################################
# Function to generate p, probability of getting an event in each group #
#########################################################################
pct_gen = function(w,mu,K,theta,tau2, M)
{
  epsi1 = matrix(rnorm(K*M,0,sqrt(0.1)), M, K)
  epsi2 = matrix(rnorm(K*M,0,sqrt(tau2)), M, K)
  pc = exp(mu+epsi1-w*(theta+epsi2))/(1+exp(mu+epsi1-w*(theta+epsi2)))
  pt = exp(mu+epsi1+(1-w)*(theta+epsi2))/(1+exp(mu+epsi1+(1-w)*(theta+epsi2)))
  return(list(pc,pt))
}
##########################################################
# Function to generate x, number of events in each group #
##########################################################
xct_gen = function(nct, pct, K)
{
  xc = matrix(NA, M, K)
  xt = matrix(NA, M, K)
  for(i in 1:K){
    xc[,i] <- sapply(1:M, function(j) rbinom(1, nct[[1]][j,i], pct[[1]][j,i]))
    xt[,i] <- sapply(1:M, function(j) rbinom(1, nct[[2]][j,i], pct[[2]][j,i]))
  }
  return(list(xc, xt))
}
####################################
# Estimators needed for some tests #
####################################
######
# DL #
######
DL <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  tau2.hat <- max(0, (Q-(sum(wt*s2)-sum(wt^2*s2)/sum(wt)))/(sum(wt)-sum(wt^2)/sum(wt)))
  return(tau2.hat)
}
######
# ML #
######
ML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  para.init <- c(sum(wt*theta.obs)/sum(wt), DL(xc, xt, nc, nt, K))
  f <- function(para) (K/2)*log(2*pi)+0.5*sum(log(s2+para[2]))+
    0.5*sum((theta.obs-para[1])^2/(s2+para[2]))
  tau2.hat <- optim(para.init, f, lower=c(-10, 0), upper=c(10,50), method="L-BFGS-B")$par[2]
  return(tau2.hat)
}
########
# REML #
########
REML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  para.init <- c(sum(wt*theta.obs)/sum(wt), DL(xc, xt, nc, nt, K))
  f <- function(para) (K/2)*log(2*pi)+0.5*sum(log(s2+para[2])) + 
    0.5*sum((theta.obs-para[1])^2/(s2+para[2])) + 0.5*log(sum(1/(s2+para[2])))
  tau2.hat <- optim(para.init, f, lower=c(-10, 0), upper=c(10,50), method="L-BFGS-B")$par[2]
  return(tau2.hat)
}
#########################
# Tests for homogeneity #
#########################
###########################
# Q-statistic based tests #
###########################
##################
# 1. Gart
##################
Q.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  P.value <- pchisq(Q, K-1, lower.tail = F)
  return(P.value)
}
#######################
# 2. Modified  Q-test #
#######################
MQ.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  m <- (nc*nt)/(nc+nt)
  v <- nc+nt-2
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  r <- 1-2/v
  w_ri <- r*wt
  theta.hat <- sum(w_ri*theta.obs)/sum(w_ri)
  Q <- sum(w_ri*(theta.obs-theta.hat)^2)
  P.value <- pchisq(Q, K-1, lower.tail = F)
  return(P.value)
}
##############
# 3. Jackson #
##############
Jackson.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/sqrt(s2)
  W <- diag(wt)
  w <- matrix(wt, ncol = 1)
  w_plus <- sum(wt)
  A <- W-(1/w_plus)*w %*% t(w)
  Y <- matrix(theta.obs, ncol = 1)
  Q <- t(Y) %*% A %*% Y
  Sigma <- diag(s2)
  S <- (Sigma^0.5) %*% A %*% Sigma^0.5
  lambda <- eigen(S)$values
  P.value <- farebrother(Q, lambda[-K])$Qq
  return(P.value)
}
############
# 4. Bliss #
############
Bliss <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  n.bar <- sum(nc+nt-2)/K
  QB <- (K-1)+sqrt((n.bar-4)/(n.bar-1))*((n.bar-2)*Q/n.bar-(K-1))
  P.value <- pchisq(QB, K-1, lower.tail = F)
  return(P.value)
}
#################
# 5. Bhaumik T3 #
#################
Bh.T3 <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  Var.Q0 <- 2*sum(wt*(s2+1/sum(wt)-2*wt*s2/sum(wt))^2)
  Var.lnQ <- Var.Q0*(1/K-1)^2
  tau.hat <- DL(xc, xt, nc, nt, K)
  wt.new <- 1/(s2+tau.hat)
  theta.hat.new <- sum(wt.new*theta.obs)/sum(wt.new)
  Q.new <- sum(wt.new*(theta.obs-theta.hat.new)^2)
  T3 <- (K-1)*(log(Q.new)-log(K-1))/sqrt(Var.lnQ)
  P.value <- pnorm(T3, lower.tail = F)
  return(P.value)
}
#######################
# 6. Kulinskaya 2015  #
#######################
# Codes from Kulinskaya, E., & Dollinger, M. B. (2015). An accurate test for homogeneity of odds ratios based on Cochran's Q-statistic. BMC medical research methodology, 15(1), 49.
source("LOR_moments_final.txt")
Q_gamma.test<-function(XC,XT,nC,nT,K){
  # nT, nC, XT, XC are vectors of length K (number of studies and number of events in T and C arms)
  a<-0
  b<-1-a
  Q<-Q_std<-numeric(0)
  EhatQ<-numeric(0)
  part<-numeric(0)
  par<-numeric(0)
  ######common theoretical parameters
  pT<-numeric(K)
  pC<-numeric(K)
  #################################
  N<-nC+nT
  qq<-nC/N  
  lenDel <- 1
  levelQ <-gam.pvalue3a.est<-0
  level1 <- 0
  pBD<-numeric(0)
  pBDT<-numeric(0)
  gamma.level3a.est1<-0
  miss<-0
  miss4a<-0
  X<-numeric(0)
  ########################################
  ##### Getting rid of zero margins ######
  #######################################
  X<-XC[(XT+XC>0)&(XT+XC<nT+nC)]
  nC0<-nC[(XT+XC>0)&(XT+XC<nT+nC)]
  nT0<-nT[(XT+XC>0)&(XT+XC<nT+nC)]
  XT<-XT[(XT+XC>0)&(XT+XC<nT+nC)]
  XC<-X
  K0<-length(nT0)
  #######################################
  ############   FOR A CHOSEN a ######################
  par<-paras_QLOR(XT,nT0,XC,nC0,a=a)
  par_std<-paras_QLOR_std(XT,nT0,XC,nC0,a=a,only0=TRUE) #### to use for standard Q
  pare<-numeric(0)
  m<-matrix(par,nrow=5)
  m1<-matrix(par_std,nrow=5)
  Deltah<-m[1,]
  wh<-m[5,]
  zetah<-m[2,]
  q<-m[3,]
  N<-m[4,]
  Q<-Q_LOR(par)[1]
  Q_std<-Q_LOR(par_std)[1]
  levelQ<-pchisq(Q_std,K0-1,ncp=0,lower.tail=FALSE)
  flatdat<-numeric(0)
  for (k in 1:K0)
  {flatdat<-c(flatdat,XT[k],XC[k],nT[k]-XT[k],nC[k]-XC[k])}
  flatdat<-array(flatdat,dim = c(2, 2, K0))
  # abd<-breslowday.test(flatdat)
  # BD<-abd$X2.HBD
  # BDT<-abd$X2.HBDT
  # pBD<-abd$p1 
  # pBDT<-abd$p2
  Deltabar<-Q_LOR(par)[2]  
  #######################################
  #####  RESTRICTIONS ON pT0, pC0????????
  #######################################
  wtm<-((2+2*cosh(b*Deltabar+zetah))/(nT0+1)+(2+2*cosh(zetah-a*Deltabar))/(nC0+1))^(-1)
  pT0<-1/(1+exp(-b*Deltabar-zetah))
  pC0<-1/(1+exp(a*Deltabar-zetah))
  #################################################
  for (k in 1:K0)
    pare<-c(pare,Deltabar,zetah[k],q[k],N[k],wtm[k]) #estimated null parametrs
  cmo<-cmomLOR_th(nT0,nC0,pT0,pC0,a=a) #### replaced cmomLOR to have exact moments
  EhatQ<-EQ(pare,cmo,nT0,nC0,a=a)
  ##### new variances and parameters calculated from models fitted through simulations
  EhatQa=0.687*EhatQ+(K0-1)*(.313) #################  new "reduced" mean as in MBD 19/11/2013
  VAR3hatQa<- -12.174*EhatQa +4.7414*(K0-1)+9.4236*EhatQa^2/(K0-1) # quadratic fit
  beta3a.est<-VAR3hatQa/EhatQa       ### adding instead Gamma appr based on corrected EQa
  alfa3a.est<-EhatQa/beta3a.est
  if(VAR3hatQa>0) gam.pvalue3a.est<-1-pgamma(Q,alfa3a.est,scale=beta3a.est) else {gam.pvalue3a.est<-(-1)
  miss4a=miss4a+1}
  ret<-c(round(gam.pvalue3a.est,3))
  #names(ret)<-c("K","Q","p_Q", "E_th(Q)", "E(Q)", "Var(Q)","alpha","beta", "p_gam", "miss", "BD stat", "p_BD", "BDT stat","p_BDT")
  return(ret)
}
############
# LR tests #
############
#############
# 7. LRT_ML #
#############
LRT.ML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt0 <- 1/s2
  theta.hat0 <- sum(wt0*theta.obs)/sum(wt0)
  tau2.hat.ML <- ML(xc, xt, nc, nt, K)
  wt.ML <- 1/(s2+tau2.hat.ML)
  theta.hat.ML <- sum(wt.ML*theta.obs)/sum(wt.ML)
  L0 <- -0.5*sum(log(s2))-0.5*sum((theta.obs-theta.hat0)^2/s2)
  L1 <- -0.5*sum(log(s2+tau2.hat.ML))-0.5*sum((theta.obs-theta.hat.ML)^2/(s2+tau2.hat.ML))
  tt <- -2*(L0-L1)
  P.value <- 0.5*pchisq(tt, 1, lower.tail = F)
  return(P.value)
}
###############
# 8. LRT_REML #
###############
LRT.REML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt0 <- 1/s2
  theta.hat0 <- sum(wt0*theta.obs)/sum(wt0)
  tau2.hat.REML <- REML(xc, xt, nc, nt, K)
  tau2.hat.ML <- ML(xc, xt, nc, nt, K)
  wt.ML <- 1/(s2+tau2.hat.ML)
  theta.hat.ML <- sum(wt.ML*theta.obs)/sum(wt.ML)
  wt.REML <- 1/(s2+tau2.hat.REML)
  L0 <- -0.5*sum(log(s2))-0.5*sum((theta.obs-theta.hat0)^2/s2)-0.5*log(sum(1/s2))
  L1 <- -0.5*sum(log(1/wt.REML))-0.5*sum(wt.REML*(theta.obs-theta.hat.ML)^2)-0.5*log(sum(wt.REML))
  tt <- -2*(L0-L1)
  P.value <- 0.5*pchisq(tt, 1, lower.tail = F)
  return(P.value)
}
###########################
# 9. Unconditional LR FeL # 
###########################
LRT.U.FE <- function(xc, xt, nc, nt, K){
  UM.RS <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="UM.FS", method="FE", nAGQ=1, to="all", drop00=FALSE)
  P.value <- UM.RS$QEp.LRT
  return(P.value)
}
############################
# 10. Unconditional LR ReL # 
############################
LRT.U.RE <- function(xc, xt, nc, nt, K){
  UM.RS <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="UM.RS", nAGQ=1, to="all", drop00=FALSE)
  P.value <- UM.RS$QEp.LRT
  return(P.value)
}
##########################
# 11. Conditional LR FeH #
##########################
LRT.C <- function(xc, xt, nc, nt, K){
  CM.RS <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="CM.AL", method = "FE", to="all", drop00=FALSE)
  P.value <- CM.RS$QEp.LRT
  return(P.value)
}
###############
# Score tests #
###############
################
# 12. Score_ML #
################
Score_ML.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  SS <- sqrt(0.5)*sum(wt^2*((theta.obs-theta.hat)^2-s2))/sqrt(sum(wt^2))
  P.value <- pnorm(SS, lower.tail = F)
  return(P.value)
}
##################
# 13. Score_REML #
##################
Score_REML.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  SS <- sqrt(0.5)*sum(wt^2*((theta.obs-theta.hat)^2-s2+1/sum(wt)))/
    sqrt(sum(wt^2)-2*sum(wt^3)/sum(wt)+sum(wt^2)^2/sum(wt)^2)
  P.value <- pnorm(SS, lower.tail = F)
  return(P.value)
}
###############################
# 14. Unconditional Score FeL #
###############################
US.FE.test <- function(xc, xt, nc, nt, K){
  Fit <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="UM.FS", method="FE", to="all", drop00=FALSE, nAGQ=1)
  Est.OR <- as.numeric(exp(Fit$beta))
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  X <- xc+xt
  aa <- Est.OR-1
  bb <- -(X+nt)*Est.OR-(nc-X)
  cc <- X*nt*Est.OR
  Delta <- bb^2-4*aa*cc
  E <- (-bb-sqrt(Delta))/(2*aa)
  VarE <- (1/E+1/(X-E)+1/(nt-E)+1/(nc-X+E))^(-1)
  TS <- sum((xt-E)^2/VarE)
  P.value <- pchisq(TS, K-1, lower.tail = F)
  return(P.value)
}
###############################
# 15. Unconditional Score ReL #
###############################
US.RE.test <- function(xc, xt, nc, nt, K){
  Fit <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="UM.RS", to="all", drop00=FALSE, nAGQ=1)
  Est.OR <- as.numeric(exp(Fit$beta))
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  X <- xc+xt
  aa <- Est.OR-1
  bb <- -(X+nt)*Est.OR-(nc-X)
  cc <- X*nt*Est.OR
  Delta <- bb^2-4*aa*cc
  E <- (-bb-sqrt(Delta))/(2*aa)
  VarE <- (1/E+1/(X-E)+1/(nt-E)+1/(nc-X+E))^(-1)
  TS <- sum((xt-E)^2/VarE)
  P.value <- pchisq(TS, K-1, lower.tail = F)
  return(P.value)
}
#############################
# 16. Conditional Score FeH #
#############################
CS.FE <- function(xc, xt, nc, nt, K){
  Fit <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="CM.AL", method = "FE", to="all", drop00=FALSE)
  Est.OR <- as.numeric(exp(Fit$beta))
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  X <- xc+xt
  ##############
  # Asymptotic
  ###############
  # aa <- Est.OR-1
  # bb <- -(X+nt)*Est.OR-(nc-X)
  # cc <- X*nt*Est.OR
  # Delta <- bb^2-4*aa*cc
  # E <- (-bb-sqrt(Delta))/(2*aa)
  # VarE <- (1/E+1/(X-E)+1/(nt-E)+1/(nc-X+E))^(-1)
  # TS <- sum((xt-E)^2/VarE)
  ##########
  # Exact
  ##########
  Hyper <- sapply(1:K, function(k) dnoncenhypergeom(NA,nt[k],nc[k],X[k],Est.OR))
  E <- sapply(1:K, function(k) sum(Hyper[[k]][,1]*Hyper[[k]][,2]))
  VarE <- sapply(1:K, function(k) sum((Hyper[[k]][,1]-E[k])^2*Hyper[[k]][,2]))
  TS <- sum((xt-E)^2/VarE)
  P.value <- pchisq(TS, K-1, lower.tail = F)
  return(P.value)
}
#####################################################
# 17. Normal Approximation to Conditional Score FeH #   
#####################################################
CS.FE.N <- function(xc, xt, nc, nt, K){
  Fit <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="CM.AL", method = "FE", to="all", drop00=FALSE)
  Est.OR <- as.numeric(exp(Fit$beta))
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  X <- xc+xt
  ##########
  # Exact
  ##########
  Hyper <- lapply(1:K, function(k) dnoncenhypergeom(NA,nt[k],nc[k],X[k],Est.OR))
  E <- sapply(1:K, function(k) sum(Hyper[[k]][,1]*Hyper[[k]][,2]))
  VarE <- sapply(1:K, function(k) sum((Hyper[[k]][,1]-E[k])^2*Hyper[[k]][,2]))
  TS <- sum((xt-E)^2/VarE)
  E2 <- sapply(1:K, function(k) sum(Hyper[[k]][,1]^2*Hyper[[k]][,2]))
  E3 <- sapply(1:K, function(k) sum(Hyper[[k]][,1]^3*Hyper[[k]][,2]))
  E4 <- sapply(1:K, function(k) sum(Hyper[[k]][,1]^4*Hyper[[k]][,2]))
  E4.Central <- sapply(1:K, function(k) sum((Hyper[[k]][,1]-E[k])^4*Hyper[[k]][,2]))
  B <- sum((E3-3*E2*E+2*E^3)/(Est.OR*VarE))
  VarTS <- sum(E4.Central/VarE^2)-K
  V <- Est.OR^2/sum(VarE)
  T4 <- (TS-K)/(VarTS-B^2*V)^(1/2)
  P.value <- pnorm(T4, lower.tail = F)
  return(P.value)
}
#############################
# 18. Conditional Score ReH #
#############################
CS.RE <- function(xc, xt, nc, nt, K){
  Fit <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="CM.AL", method = "ML", to="all", drop00=FALSE)
  Est.OR <- as.numeric(exp(Fit$beta))
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  X <- xc+xt
  ##########
  # Exact
  ##########
  Hyper <- lapply(1:K, function(k) dnoncenhypergeom(NA,nt[k],nc[k],X[k],Est.OR))
  E <- sapply(1:K, function(k) sum(Hyper[[k]][,1]*Hyper[[k]][,2]))
  VarE <- sapply(1:K, function(k) sum((Hyper[[k]][,1]-E[k])^2*Hyper[[k]][,2]))
  E3.Central <- sapply(1:K, function(k) sum((Hyper[[k]][,1]-E[k])^3*Hyper[[k]][,2]))
  S <- sum((xt-E)^2-VarE)
  I11 <- sum(sapply(1:K, function(k) sum(((Hyper[[k]][,1]-E[k])^2-VarE[k])^2*Hyper[[k]][,2])))
  I12 <- sum(E3.Central)
  I22 <- sum(VarE)
  T5 <- S/sqrt(I11-I12^2/I22)
  P.value <- pnorm(T5, lower.tail = F)
  return(P.value)
}
##########
# 19. BD #
##########
BD <- function(xc, xt, nc, nt, K){
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  Est.OR <- sum(xt*(nc-xc)/(nc+nt))/sum(xc*(nt-xt)/(nc+nt))
  X <- xc+xt
  aa <- Est.OR-1
  bb <- -(X+nt)*Est.OR-(nc-X)
  cc <- X*nt*Est.OR
  Delta <- bb^2-4*aa*cc
  E <- (-bb-sqrt(Delta))/(2*aa)
  VarE <- (1/E+1/(X-E)+1/(nt-E)+1/(nc-X+E))^(-1)
  TS <- sum((xt-E)^2/VarE)
  P.value <- pchisq(TS, K-1, lower.tail = F)
  return(P.value)
}
###########
# 20. MBD #
###########
MBD <- function(xc, xt, nc, nt, K){
  xc <- xc+0.5
  xt <- xt+0.5
  nc <- nc+1
  nt <- nt+1
  Est.OR <- sum(xt*(nc-xc)/(nc+nt))/sum(xc*(nt-xt)/(nc+nt))
  X <- xc+xt
  aa <- Est.OR-1
  bb <- -(X+nt)*Est.OR-(nc-X)
  cc <- X*nt*Est.OR
  Delta <- bb^2-4*aa*cc
  E <- (-bb-sqrt(Delta))/(2*aa)
  VarE <- (1/E+1/(X-E)+1/(nt-E)+1/(nc-X+E))^(-1)
  TS <- sum((xt-E)^2/VarE)-(sum(xt)-sum(E))/sum(VarE)
  P.value <- pchisq(TS, K-1, lower.tail = F)
  return(P.value)
}
##############
# Wald tests #
##############
###############
# 21. Wald_ML #
###############
Wald.ML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  tau2.hat.ML <- ML(xc, xt, nc, nt, K)
  wt.ML <- 1/(s2+tau2.hat.ML)
  Var.ML <- 2/sum(wt.ML^2)
  tt <- tau2.hat.ML/sqrt(Var.ML)
  P.value <- pnorm(tt, lower.tail = F)
  return(P.value)
}
#################
# 22. Wald_REML #
#################
Wald.REML <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  tau2.hat.REML <- REML(xc, xt, nc, nt, K)
  wt.REML <- 1/(s2+tau2.hat.REML)
  Var.REML <- 2/(sum(wt.REML^2)-2*sum(wt.REML^3)/sum(wt.REML)+sum(wt.REML^2)^2/sum(wt.REML)^2)
  tt <- tau2.hat.REML/sqrt(Var.REML)
  P.value <- pnorm(tt, lower.tail = F)
  return(P.value)
}
###############
# 23. Wald_FE #
###############
# Do not work in most of the rare binary events settings
# Wald.FE <- function(xc, xt, nc, nt, K){
#   # UM.RS <- rma.glmm(xt, nt-xt, xc, nc-xc, measure = "OR", model="UM.FS", method="FE", nAGQ=1, to="all", drop00=FALSE)
#   # P.value <- UM.RS$QEp.Wld
#   Y <- as.factor(unlist(lapply(1:K, function(k) c(rep(1, xt[k]), rep(0,nt[k]-xt[k]), rep(1, xc[k]), rep(0,nc[k]-xc[k])))))
#   Cluster <- as.factor(unlist(lapply(1:K, function(k) rep(k,nt[k]+nc[k]))))
#   Treat <- as.factor(unlist(lapply(1:K, function(k) c(rep(1,nt[k]), rep(0,nc[k])))))
#   Fit <- glm(Y~1+Treat+Cluster+Treat*Cluster, family=binomial)
#   P.value <- as.numeric(wald.test(vcov(Fit),coef(Fit), Terms=(2+K):(2*K))$result$chi2[3])
#   return(P.value)
# }
###############
# Other tests #
###############
############
# 24. Peto #
############
Peto <- function(xc, xt, nc, nt, K){
  Fit <- rma.peto(xt, nt-xt, xc, nc-xc, to="all", drop00=FALSE)
  P.value <- Fit$QEp
  return(P.value)
}
###############
# 25. Z^2_WLS #
###############
Z2_WLS.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  Z2_WLS <- (Q-(K-1))^2/(2*(K-1))
  P.value <- pf(Z2_WLS, 1, K-1, lower.tail = F)
  return(P.value)
}
################
# 26. Z^2_WLSR #
################
Z2_WLSR.test <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(xt+0.5)+1/(nt-xt+0.5)+1/(xc+0.5)+1/(nc-xc+0.5)
  wt <- 1/s2
  theta.hat <- sum(wt*theta.obs)/sum(wt)
  Q <- sum(wt*(theta.obs-theta.hat)^2)
  Z2_WLSR <- (Q-(K-1))^2/sum(((theta.obs-theta.hat)^2*wt-1)^2)
  P.value <- pf(Z2_WLSR, 1, K-1, lower.tail = F)
  return(P.value)
}
#############
# 27. Z^2_K #
#############
Z2_K.test <- function(xc, xt, nc, nt, K){
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(xt+0.5)+1/(nt-xt+0.5)+1/(xc+0.5)+1/(nc-xc+0.5)
  theta.hat <- sum(theta.obs/(1/nc+1/nt))/sum(1/(1/nc+1/nt))
  V1 <- ((theta.obs-theta.hat)^2-s2)/(1/nc+1/nt)^2
  V2 <- ((theta.obs-theta.hat)^2-s2)^2/(1/nc+1/nt)^4
  Z2_K <- sum(V1)^2/sum(V2)
  P.value <- pf(Z2_K, 1, K-1, lower.tail = F)
  return(P.value)
}
##################
# 28. Bhaumik T4 #
##################
Bh.T4 <- function(xc, xt, nc, nt, K){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  theta.hat.sa <- mean(theta.obs)
  y <- (theta.obs-theta.hat.sa)^2
  B <- ((K-2)/K)*s2 + sum(s2)/K^2
  T4 <- sum(y-B)/sqrt(sum(2*B))
  P.value <- pnorm(T4, lower.tail = F)
  return(P.value)
}
#######################################
# 29. Parametric bootstrap Bhaumik T4 #
#######################################
Bh.T4.Boot <- function(xc, xt, nc, nt, K, B.time=200){
  pc <- (xc+0.5)/(nc+1)
  pt <- (xt+0.5)/(nt+1)
  s2 <- 1/(nt*pt*(1-pt))+1/(nc*pc*(1-pc))
  theta.obs <- log((xt+0.5)/(nt-xt+0.5)) - log((xc+0.5)/(nc-xc+0.5))
  mu.obs <- log((xc+0.5)/(nc-xc+0.5))
  theta.hat.sa <- mean(theta.obs)
  pc.hat <- pc
  pt.hat <- exp(mu.obs+theta.hat.sa)/(1+exp(mu.obs+theta.hat.sa))
  y <- (theta.obs-theta.hat.sa)^2
  B <- ((K-2)/K)*s2 + sum(s2)/K^2
  T4 <- sum(y-B)/sqrt(sum(2*B))
  T4.b <- rep(NA, B.time)
  for(i in 1:B.time){
    xc.b <- sapply(1:K, function(j) rbinom(1,nc[j],pc.hat[j]))
    xt.b <- sapply(1:K, function(j) rbinom(1,nt[j],pt.hat[j]))
    pc.b <- (xc.b+0.5)/(nc+1)
    pt.b <- (xt.b+0.5)/(nt+1)
    theta.obs.b <- log((xt.b+0.5)/(nt-xt.b+0.5)) - log((xc.b+0.5)/(nc-xc.b+0.5))
    s2.b <- 1/(nt*pt.b*(1-pt.b))+1/(nc*pc.b*(1-pc.b))
    theta.hat.sa.b <- mean(theta.obs.b)
    y.b <- (theta.obs.b-theta.hat.sa.b)^2
    B.b <- ((K-2)/K)*s2.b + sum(s2.b)/K^2
    T4.b[i] <- sum(y.b-B.b)/sqrt(sum(2*B.b))
  }
  P.value <- sum(T4.b>T4)/B.time
  return(P.value)
}
##################################
# 30. Nonparametric Bootstrap DL #
##################################
DLb.test <- function(xc, xt, nc, nt, K, B=200){
  id <- matrix(sample(1:K, K*B, replace = T), B, K)
  tau2.hats <- rep(NA, B)
  for(i in 1:B){
    xc.b <- xc[id[i,]]
    xt.b <- xt[id[i,]]
    nc.b <- nc[id[i,]]
    nt.b <- nt[id[i,]]
    pc <- (xc.b+0.5)/(nc.b+1)
    pt <- (xt.b+0.5)/(nt.b+1)
    theta.obs <- log((xt.b+0.5)/(nt.b-xt.b+0.5)) - log((xc.b+0.5)/(nc.b-xc.b+0.5))
    s2 <- 1/(nt.b*pt*(1-pt))+1/(nc.b*pc*(1-pc))
    wt <- 1/s2
    theta.hat <- sum(wt*theta.obs)/sum(wt)
    Q <- sum(wt*(theta.obs-theta.hat)^2)
    tau2.hats[i] <- max(0, (Q-(sum(wt*s2)-sum(wt^2*s2)/sum(wt)))/(sum(wt)-sum(wt^2)/sum(wt)))
  }
  P.value <- sum(tau2.hats==0)/B
  return(P.value)
}

# Put all compared methods into a list
Methods <- list(Q.test, MQ.test, Jackson.test, Bliss, Bh.T3, Q_gamma.test, LRT.ML, LRT.REML, LRT.U.FE, LRT.U.RE, LRT.C, 
                Score_ML.test, Score_REML.test, US.FE.test, US.RE.test, CS.FE, CS.FE.N, CS.RE, BD, MBD, Wald.ML, 
                Wald.REML, Peto, Z2_WLS.test, Z2_WLSR.test, Z2_K.test, Bh.T4, Bh.T4.Boot, DLb.test)
No.Methods <- length(Methods)

# Function to obtain p-values under null hypothesis
Emp_Dist <- function(low,up,R,w,M,theta,K,mu,tau2=0){
  N <- nct_gen(low, up, R, K, M)
  P <- pct_gen(w, mu, K, theta, tau2, M)
  X <- xct_gen(N, P, K)
  P.values <- matrix(NA, M, No.Methods)
  for(i in 1:No.Methods){
    P.values[,i] <- sapply(1:M, function(j) Methods[[i]](X[[1]][j,], X[[2]][j,], N[[1]][j,], N[[2]][j,], K))
  }
  return(P.values)
}

##############################
# Example for implementation #
##############################
################################################
# R=1, K=20, mu=-5, theta=0, small sample size #
################################################
R <- 1
low <- 20
up <- 1000
K <- 20
mu <- -5
theta <- 0
M <- 1000
ED1 <- Emp_Dist(low,up,R,0, M, theta, K, mu, tau2=0)
ED2 <- Emp_Dist(low,up,R,0.5, M, theta, K, mu, tau2=0)
ED3 <- Emp_Dist(low,up,R,1, M, theta, K, mu, tau2=0)

# Critical value to maintain type I error rate at 0.05
Emp_Cut1 <- apply(ED1, 2, quantile, 0.05)
Emp_Cut2 <- apply(ED2, 2, quantile, 0.05)
Emp_Cut3 <- apply(ED3, 2, quantile, 0.05)


tau2 <- seq(0.1, 1, 0.1)
# Function to obtain adjusted power
Comparison.Test <- function(low,up,R,w,M,theta,K,mu,tau2, Emp_Cut){
  N <- nct_gen(low, up, R, K, M)
  P <- pct_gen(w, mu, K, theta, tau2, M)
  X <- xct_gen(N, P, K)
  P.values <- matrix(NA, M, No.Methods)
  Adj.Power <- rep(NA, No.Methods)
  for(i in 1:No.Methods){
    P.values[,i] <- sapply(1:M, function(j) Methods[[i]](X[[1]][j,], X[[2]][j,], N[[1]][j,], N[[2]][j,], K))
    # Adjusted power is calculated by using the empirical critical value to maintain Type I error rate exactly 0.05
    Adj.Power[i] <- ecdf(P.values[,i])(Emp_Cut[i])  
  }
  return(Adj.Power)
}

Estimate1 = lapply(tau2, function(x) Comparison.Test(low,up,R,0, M, theta, K, mu, x, Emp_Cut1))
Estimate2 = lapply(tau2, function(x) Comparison.Test(low,up,R,0.5, M, theta, K, mu, x, Emp_Cut2))
Estimate3 = lapply(tau2, function(x) Comparison.Test(low,up,R,1, M, theta, K, mu, x, Emp_Cut3))
