
loglikLOOCVVAR1sim <- function (lambdas, Y, penalty="ridge", ...){
    if (as.character(class(Y)) != "array") {
        stop("Input (Y) is of wrong class.")
    }
    if (length(dim(Y)) != 3) {
        stop("Input (Y) is of wrong dimensions: either covariate, time or sample dimension is missing.")
    }
    if (as.character(class(lambdas)) != "numeric") {
        stop("Input (lambdas) is of wrong class.")
    }
    if (length(lambdas) != 2) {
        stop("Input (lambdas) is of wrong length.")
    }
    if (any(is.na(lambdas))) {
        stop("Input (lambdas) is not a vector of non-negative numbers.")
    }
    if (any(lambdas < 0)) {
        stop("Input (lambdas) is not a vector of non-negative numbers.")
    }
    loglik <- 0
    
    # ridge leave K-fold-out cross-validated log-likelihood
    if(penalty=="ridge"){
              for (k in 1:dim(Y)[3]) {
                        # Ridge ML eximtion of the VAR(1) model
                        VAR1hat <- ridgeVAR1(Y[,,-k], lambdas[1], lambdas[2], ...)
                        for (l in 2:dim(Y)[2]){
                                  res <- Y[, l, k] - VAR1hat$A %*%Y[, l-1, k]
                                  loglik <- loglik - t(res) %*% VAR1hat$P %*% res/2 + log(det(VAR1hat$P))/2
                        }
              }
              return(-loglik)
    }
    # SCAD leave K-fold-out cross-validated log-likelihood
    if(penalty=="scad"){
              for (k in 1:dim(Y)[3]) {
                        # convert the data without k-th cell line
                        Ylaso <- array2longitudinal(Y[,,-k])
                        Ylaso1 <- as.longitudinal(Ylaso,repeats=dim(Y)[3]-1)   
                        # SCAD estimation of VAR(1) model                   
                        VAR1hat <- sparse.tscgm(data=Ylaso1, lam1=lambdas[1], lam2=lambdas[2], model="ar1", optimality=NULL, control = list(maxit.out=5, maxit.in=10))
                        Omega <-VAR1hat$theta

                        for (l in 2:dim(Y)[2]){
                                  res <- Y[, l, k] - t(VAR1hat$gamma) %*%Y[, l-1, k]
                                  loglik <- loglik - t(res) %*% Omega %*% res/2 + log(det(Omega))/2
                        }
              }
              return(-loglik)
    }
}
