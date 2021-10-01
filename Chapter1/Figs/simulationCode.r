\newpage


\begin{lstlisting}[language=R]
loglikLOOCVVAR1sim <- function (lambdas, Y, penalty="ridge", ...){
    loglik <- 0
    if(penalty=="ridge"){
       for (k in 1:dim(Y)[3]) {
          VAR1hat <- ridgeVAR1(Y[,,-k], lambdas[1], lambdas[2], ...)
           for (l in 2:dim(Y)[2]){
              res <- Y[, l, k] - VAR1hat$A %*%Y[, l-1, k]
              loglik <- loglik - t(res) %*% VAR1hat$P %*% res/2 + 
	                   log(det(VAR1hat$P))/2}}
              return(-loglik)}
    if(penalty=="scad"){
       for (k in 1:dim(Y)[3]) {
         Ylaso <- array2longitudinal(Y[,,-k])
         Ylaso1 <- as.longitudinal(Ylaso,repeats=dim(Y)[3]-1)                                                          
         VAR1hat <- sparse.tscgm(data=Ylaso1, lam1=lambdas[1], lam2=lambdas[2], 
				optimality=NULL)
         for (l in 2:dim(Y)[2]){
            res <- Y[, l, k] - t(VAR1hat$gamma) %*%Y[, l-1, k]
            loglik <- loglik - t(res) %*% VAR1hat$theta %*% res/2 +
			log(det(VAR1hat$theta))/2}}
   return(-loglik)}
}


simDataGen< function(p,n,T,data="sim",topologyA="clique",topologyP="banded"){
  if(data=="sim"){
    trueA <- createA(p,topology=topologyA,nonzeroA=0.3, nCliques=8, 
	nHubs=p%/%10, percZeros=0.95, stationary=TRUE)
    trueP<-createS(n,p,topology=topologyP,nonzero=0.5,banded.n =3,precision = T)}
  if(data=="real"){
      data(hpvP53)
      GEdat <- exprs(hpvP53)
      Yrel <- array(NA, dim = c(p, 8, 4))
      Yrel[,,1] <- GEdat[1:p,1:8]
      Yrel[,,2] <- GEdat[1:p,9:16]
      Yrel[,,3] <- GEdat[1:p,17:24]
      Yrel[,,4] <- GEdat[1:p,25:32]
      VAR1hat <- ridgeVAR1(Yrel, 1, 0.1)
      trueA <-VAR1hat$A   # trueA1 for sparse comparison
      trueP <-VAR1hat$P}
   if(data=="realSpar"){
      trueA <- matrix(0, p, p)
      trueA[sparsifyVAR1(A=trueA1,SigmaE=symm(solve(trueP)),threshold="localFDR",
	FDRcut=0.8)$nonzeros] <- 0.5
      trueP=sparsify(trueP, threshold = "localFDR",FDRcut = 0.8)[[2]]
      trueP[trueP == 1] <- 0.5}
   return(list(A=trueA,P=trueP))
}
\end{lstlisting}

\newpage

\begin{lstlisting}[language=R]
rm(list=ls())
set.seed(321)
# load libraries
library(longitudinal)
library(Biobase)
library(SparseTSCGM)
library(ragt2ridges)
setwd("N:/Cervical Cancer Data/Patways/Review")
source("loglikLOOCVVAR1sim.R")
source("simDataGen.R")

p=25 # numbers of variables (genes)
n=5  # number of samples (cell lines)
T=10 # number of time points

#  Generate true paramterer matrix and precision matrix
sim <- simDataGen(p=p,n=n,T=T,data="sim",topologyA="clique", topologyP="banded")
trueA <- sim$A;trueP <- sim$P

#  Generate and center data
Y <- dataVAR1(n, T, trueA, solve(trueP))
Yridge <- centerVAR1data(Y)

#  ragt2ridges estimator
#  using optimal penalities estimate A and Sigma
optL1 <- optPenaltyVAR1(Yridge, rep(10^(-10), 2), rep(30, 2), c(1,0.1), 
		optimizer="nlm")
rdgEst1 <- ridgeVAR1(Y=Yridge,lambdaA=optL1[1],lambdaP=optL1[2])
ridgeA1 <- rdgEst1$A; ridgeP1 <- rdgEst1$P

optL2 <- nlminb(c(optL1[1],optL1[2]), loglikLOOCVVAR1sim, gradient=NULL,
		 lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="ridge",
		 control=list(rel.tol=0.01))$par
rdgEst2 <- ridgeVAR1(Y=Yridge,lambdaA=optL2[1],lambdaP=optL2[2])
ridgeA2 <- rdgEst2$A; ridgeP2 <- rdgEst2$P

#  convert array data to longitudinal
Ylaso <- array2longitudinal(Yridge)
Ylaso1 <- as.longitudinal(Ylaso,repeats=n)

#  SparseTSCGM estimator of A and Sigma with optimal lambda
res.tscgm3 <- sparse.tscgm(data=Ylaso1, lam1=NULL, lam2=NULL, model="ar1",
		 optimality="bic")
lasoA3 <- t(res.tscgm3$gamma); lasoP3 <- res.tscgm3$theta

#  estimate A and sigma using CV penality paramters
optL4 <- nlminb(c(res.tscgm3$lam1.opt,res.tscgm3$lam2.opt), loglikLOOCVVAR1sim,
	 gradient=NULL, lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="scad", 
	 control=list(rel.tol=0.01))$par
res.tscgm4 <- sparse.tscgm(data=Ylaso1, lam1=optL4[1], lam2=optL4[2], 
		optimality=NULL)
lasoA4 <- t(res.tscgm4$gamma); lasoP4 <- res.tscgm4$theta

#   Forbenius loss calculation for A and Sigma
lossPr <- sqrt(sum((as.numeric(ridgeP1)-as.numeric(trueP))^2))
lossPl <- sqrt(sum((as.numeric(lasoP3)-as.numeric(trueP))^2))
lossAr <- sqrt(sum((as.numeric(ridgeA1)-as.numeric(trueA))^2))
lossAl <- sqrt(sum((as.numeric(lasoA3)-as.numeric(trueA))^2))
cvlossPr <- sqrt(sum((as.numeric(ridgeP2)-as.numeric(trueP))^2))
cvlossPl <- sqrt(sum((as.numeric(lasoP4)-as.numeric(trueP))^2))
cvlossAr <- sqrt(sum((as.numeric(ridgeA2)-as.numeric(trueA))^2))
cvlossAl <- sqrt(sum((as.numeric(lasoA4)-as.numeric(trueA))^2))
\end{lstlisting}

\newpage
\begin{lstlisting}[language=R]
rm(list=ls())
set.seed(321)

# load libraries
library(longitudinal)
library(Biobase)
library(SparseTSCGM)
library(ragt2ridges)
library(ROCR)
setwd("N:/Cervical Cancer Data/Patways/Review")
source("loglikLOOCVVAR1sim.R")
source("simDataGen.R")

p=25 # numbers of variables (genes)
n=15  # number of samples (cell lines)
T=20 # number of time points

#  Generate true paramterer matrix and precision matrix        
sim <- simDataGen(p=p,n=n,T=T,data="sim",topologyA="clique", topologyP="banded")
trueA1 <- sim$A;trueP <- sim$P

trueA <-trueA1
trueA[trueA!=0] <- 1

#  Generate and center simulated data
Y <- dataVAR1(n, T, trueA1, solve(trueP))
Yridge <- centerVAR1data(Y)

#  SparseTSCGM package
Ylaso <- array2longitudinal(Yridge)
Ylaso1 <- as.longitudinal(Ylaso,repeats=n)

SCADa <- sparse.tscgm(data=Ylaso1, lam1=NULL, lam2=NULL, model="ar1",
		 optimality="bic")
lasoAa <- t(SCADa$gamma)
#  Calculate true and false postive rate
PRDa <- prediction(as.numeric(lasoAa),as.numeric(trueA))
TPRlA <- PRDa@tp[[1]][2]/max(PRDa@tp[[1]])
FPRlA <- PRDa@fp[[1]][2]/max(PRDa@fp[[1]])

optLlb <- nlminb(c(SCADa$lam1.opt,SCADa$lam2.opt), loglikLOOCVVAR1sim, 
	gradient=NULL, lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="scad", 
	control=list(rel.tol=0.01))$par
lasoAb <- t(sparse.tscgm(data=Ylaso1, lam1=optLlb[1], lam2=optLlb[2],
	 optimality=NULL)$gamma)
PRDb <- prediction(as.numeric(lasoAb),as.numeric(trueA))
TPRlB <- PRDb@tp[[1]][2]/max(PRDb@tp[[1]])
FPRlB <- PRDb@fp[[1]][2]/max(PRDb@fp[[1]])

#  ragt2ridges package
optLrA <- optPenaltyVAR1(Yridge, rep(10^(-10), 2), rep(30, 2), c(1,0.1), optimizer="nlm")
rdgEstA <- ridgeVAR1(Y=Yridge,lambdaA=optLrA[1],lambdaP=optLrA[2])
##################       a      ##################################
ridgeAa <- matrix(0, p, p)
ridgeAa[sparsifyVAR1(A=rdgEstA$A,SigmaE=symm(solve(rdgEstA$P)),threshold="top",
		top=as.numeric(sum(lasoAa == 1)))$nonzeros] <- 1
#  Calculate false positive and true postive rate
predA <- prediction(as.numeric(ridgeAa),as.numeric(trueA))
TPRrA <- predA@tp[[1]][2]/max(predA@tp[[1]])
FPRrA <- predA@fp[[1]][2]/max(predA@fp[[1]])
##################       b      ##################################
ridgeAb <- matrix(0, p, p)
ridgeAb[sparsifyVAR1(A=rdgEstA$A, SigmaE=symm(solve(rdgEstA$P)),threshold="localFDR",
		FDRcut=0.8)$nonzeros] <- 1
#  Calculate false positive and true postive rate
predB <- prediction(as.numeric(ridgeAb),as.numeric(trueA))
TPRrB <- predB@tp[[1]][2]/max(predB@tp[[1]])
FPRrB <- predB@fp[[1]][2]/max(predB@fp[[1]])
##################       c      ##################################
optLrC <- nlminb(c(optLrA[1],optLrA[2]), loglikLOOCVVAR1sim, gradient=NULL,
	 lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="ridge", control=list(rel.tol=0.01))$par
rdgEstC <- ridgeVAR1(Y=Yridge,lambdaA=optLrC[1],lambdaP=optLrC[2])
ridgeAc <- matrix(0, p, p)
ridgeAc[sparsifyVAR1(A=rdgEstC$A, SigmaE=symm(solve(rdgEstC$P)),
	threshold="localFDR",FDRcut=0.8)$nonzeros] <- 1
#  Calculate false positive and true postive rate
predC <- prediction(as.numeric(ridgeAb),as.numeric(trueA))
TPRrC <- predC@tp[[1]][2]/max(predC@tp[[1]])
FPRrC <- predC@fp[[1]][2]/max(predC@fp[[1]])
\end{lstlisting}

