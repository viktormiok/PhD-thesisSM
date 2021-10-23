# set seed
set.seed(321)

# load libraries
library(longitudinal)
library(Biobase)
library(Matrix)
library(SparseTSCGM)
library(ragt2ridges)

# declare function
loglikLOOCVVAR1sim <- function (lambdas, Y, penalty="ridge", ...){
	loglik <- 0
    
	# ridge leave K-fold-out cross-validated log-likelihood
	if(penalty=="ridge"){
		for (k in 1:dim(Y)[3]) {
			# Ridge ML eximtion of the VAR(1) model
			VAR1hat <- ridgeVAR1(Y[,,-k], lambdas[1], lambdas[2], ...)
			for (l in 2:dim(Y)[2]){
				res <- Y[, l, k] - VAR1hat$A %*%Y[, l-1, k]
				loglik <- loglik - t(res) %*% VAR1hat$P %*% res/2 + 
						log(det(VAR1hat$P))/2
			}
		}
	}

	# SCAD leave K-fold-out cross-validated log-likelihood
	if(penalty=="scad"){
		for (k in 1:dim(Y)[3]) {
			# convert the data without k-th cell line
			Yscad <- array2longitudinal(Y[,,-k])
			Yscad1 <- as.longitudinal(Yscad,repeats=dim(Y)[3]-1)   

			# SCAD estimation of VAR(1) model                   
			VAR1hat <- sparse.tscgm(data=Yscad1, lam1=lambdas[1], 
					lam2=lambdas[2], model="ar1", optimality=NULL, 
					control=list(maxit.out=5, maxit.in=10))
			Omega <- VAR1hat$theta
			for (l in 2:dim(Y)[2]){
				res <- Y[, l, k] - t(VAR1hat$gamma) %*%Y[, l-1, k]
				loglik <- loglik - t(res) %*% Omega %*% res/2 + 
						log(det(Omega))/2
			}
		}
	}

	return(-loglik)
}


# set number variates, samples, time points
p=25
n=15
T=20

# generate autoregressive coefficient matrix
trueA <- createA(p, "clique", 0.3, 8)

# generate precision matrix
trueP <- createS(n, p, "banded", nonzero=0.5, banded.n=3, precision=T)

# sample data from a VAR(1) model
Y <- dataVAR1(n, T, trueA, solve(trueP))

# covariate-wise zero centering of the data
Yridge <- centerVAR1data(Y)

# optimal lambdas for ragt2ridges
lamR<-seq(4,0.01,length.out=20)
LOOCVrPrev <- loglikLOOCVVAR1sim(c(lamR[1],lamR[1]), Yridge,
					 penalty="ridge")
LOOCVr <- numeric()
for(i in 1:length(lamR)) {
	LOOCVr <- cbind(LOOCVr,loglikLOOCVVAR1sim(c(lamR[i],lamR[i]), 
			Yridge, penalty="ridge"))
	if (LOOCVr[i]-LOOCVrPrev > 5){
		k <- which.min(LOOCVr)
		break
	} else {
		LOOCVrPrev = LOOCVr[i]
		k <- which.min(LOOCVr)
	}
}
optLr <- nlminb(c(lamR[k], lamR[k]), loglikLOOCVVAR1sim,gradient=NULL, 
		lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="ridge", 
		control=list(rel.tol=10^(-10)))$par

# Ridge ML eximtion of the VAR(1) model
ridgeEst <- ridgeVAR1(Y=Yridge, lambdaA=optLr[1], lambdaP=optLr[2])
ridgeA <- ridgeEst$A 
ridgeP <- ridgeEst$P
              
# convert a time-series array to a longitudinal object required for SparseTSCGM
Yscad <- array2longitudinal(Yridge)
Yscad1 <- as.longitudinal(Yscad,repeats=n)

# optimal lambdas for ridge method
lamS<-seq(2,0.01,length.out =20)
LOOCVs <- numeric()
LOOCVsPrev <- loglikLOOCVVAR1sim(c(lamS[1],lamS[1]), Yridge, 
					penalty="scad")
for(j in 1:length(lamS)) {
	LOOCVs <- cbind(LOOCVs,loglikLOOCVVAR1sim(c(lamS[j],lamS[j]), Yridge, 
			penalty="scad"))        
        if (LOOCVs[j]-LOOCVsPrev>5){
		l <- which.min(LOOCVs)
		break
	} else{
		LOOCVsPrev=LOOCVs[j]
		l <- which.min(LOOCVs)
	}
} 
optLl <- nlminb(c(lamS[l],lamS[l]), loglikLOOCVVAR1sim, gradient=NULL, 
		lower=c(10^(-10),10^(-10)), Y=Yridge, penalty="scad", 
		control=list(rel.tol=10^(-10)))$par

# SCAD estimation of autregressive coefficient and precision matrix
scadEst <- sparse.tscgm(data=Yscad1, lam1=optLl[1], lam2=optLl[2], optimality=NULL)
scadA <- t(scadEst$gamma)
scadP <- scadEst$theta

# Frobenius loss for A and Sigma
sqrt(sum((as.numeric(ridgeP) - as.numeric(trueP))^2))
sqrt(sum((as.numeric(scadP) - as.numeric(trueP))^2))
sqrt(sum((as.numeric(ridgeA) - as.numeric(trueA))^2))
sqrt(sum((as.numeric(scadA) - as.numeric(trueA))^2))
