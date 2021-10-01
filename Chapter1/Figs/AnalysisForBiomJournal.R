# load libraries
library(ragt2ridges)
library(Biobase)

# load and reformat data
data(hpvP53)
Y <- longitudinal2array(t(exprs(hpvP53)))

# search for optimal penalty parameters
optLambdas1 <- optPenaltyVAR1(Y, lambdaMin=c(120, 0.0024), lambdaMax=c(930, 0.0026), lambdaInit=c(250, 0.0025))

# specify grid for contour
lambdaAgrid <- seq(120, 930, length.out=20)
lambdaPgrid <- sort(c(seq(0.002557, .0026, length.out=100), seq(0.001557, .0036, length.out=20)))
LOOCVres1 <- loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, Y)

# plot contour
contour(lambdaAgrid, lambdaPgrid, LOOCVres1$llLOOCV, xlab = "lambdaA", ylab = "lambdaP", main = "cross-validated log-likelihood", nlevels=25)
points(optLambdas1[1], optLambdas1[2], pch=20, cex=2, col="red")

# fit VAR(1) model
VAR1hat <- ridgeVAR1(Y=Y, lambdaA=optLambdas1[1], lambdaP=optLambdas1[2])
Ahat <- VAR1hat$A; Phat <- VAR1hat$P
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <- colnames(Phat) <- rownames(hpvP53)
edgeHeat(Ahat, main="ridge estimate of A")
edgeHeat(Phat, main="ridge precision estimate")

# determine support for A and O
zerosA <- sparsifyVAR1(A=Ahat, SigmaE=symm(solve(Phat)), threshold="localFDR", FDRcut=0.95, statistics=F)$zeros
zerosP <- sparsify(Phat, threshold = "localFDR",FDRcut = 0.95, output="light")$zeros

# determine optimal lambda's with inferred support
optLambdas2 <- optPenaltyVAR1(Y, lambdaMin=c(10^(-5),10^(-5)), lambdaMax=c(10,0.1), lambdaInit=c(5,0.01), zerosA=zerosA, zerosP=zerosP, zerosAfit="sparse")

# determine contour
lambdaAgrid <- seq(-1, 1, length.out=20) + optLambdas2[1]
lambdaPgrid <- seq(-0.001, 0.001, length.out=20) + optLambdas2[2]
LOOCVres2 <- loglikLOOCVcontourVAR1(lambdaAgrid, lambdaPgrid, Y, zerosA=zerosA, zerosP=zerosP, zerosAfit="sparse")

# plot contour
contour(lambdaAgrid, lambdaPgrid, LOOCVres2$llLOOCV, xlab = "lambdaA", ylab = "lambdaP", main = "cross-validated log-likelihood", nlevels=25)
points(optLambdas2[1], optLambdas2[2], pch=20, cex=2, col="red")

# re-fit of the VAR(1) model including the prior knowledge
AhatNonsparse <- Ahat
PhatNonsparse <- Phat
VAR1hat <- ridgeVAR1(Y=Y, lambdaA=optLambdas2[1], lambdaP=optLambdas2[2], zerosA=zerosA, zerosP=zerosP, zerosAfit="sparse")
Ahat <- VAR1hat$A; Phat <- VAR1hat$P
rownames(Ahat) <- colnames(Ahat) <- rownames(Phat) <- colnames(Phat) <- rownames(hpvP53)
edgeHeat(Ahat, main="ridge re-estimate of A, with inferred support")
edgeHeat(Phat, main="ridge precision re-estimate, with inferred support")
edgeHeat(adjacentMat(Ahat), legend=FALSE, main="inferred support of A")
edgeHeat(adjacentMat(Phat), legend=FALSE, main="inferred support of precision matrix")

# graph of the temporal interaction among the genes
graphVAR1(Ahat, Phat, nNames=rownames(Ahat), type="TSCG",vertex.label.cex=0.5, vertex.label.font=1, vertex.size=4, vertex.label.color.T0="darkblue", 
     vertex.label.color.T1="darkblue",  vertex.frame.color="steelblue", vertex.color.T0="lightblue", 
     vertex.color.T1="lightblue", edge.width=1.5, main="")

# node statistics table
stats <- nodeStatsVAR1(Ahat, Phat, as.table=TRUE)
rownames(stats) <- fData(hpvP53)[,1]

