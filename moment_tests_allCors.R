nX <- 1
bX <- 20
rhos <- 0.4
Beta <- 0
s2 <- 1
nsims0 <- 1000
nsims <- 100
ndata <- 1000
m <- bX

library(MASS)
library(Matrix)

# Set intra-correlations of X's
rho_blocksX <- lapply(rhos, function (R) matrix(R, bX, bX) + 
                        diag(rep(1 - R, bX)))
mX <- nX * bX

#------Generating Block Covariance Matrices-----#
SigmaX <- as.matrix(bdiag(rho_blocksX))

# Setting up loop
source("moment_calcs.R")
source("varcalcs_chiSq.R")
source("mvrnormR.R")

meanMat <- matrix(0, m, m)
cr2sum <- estvar <- numeric(nsims0 * nsims)

for (sim0 in 1:nsims0) {
  
  corsMat <- matrix(0, bX, nsims)
  
  for (sim in 1:nsims) {
    
    simpos <- (sim0 - 1) * nsims + sim
    
    cat("simpos", simpos, "/", nsims * nsims0, "\n")
    
    Data <- mvrnorm(ndata, rep(0, bX), SigmaX)
    Y <- Data %*% rep(Beta, bX) + rnorm(ndata, sd = sqrt(s2))
    corsMat[ , sim] <- cor(Data, Y)
    estvar[simpos] <- varcalc1(Y, Data)
    cr2sum[simpos] <- crossprod(corsMat[ , sim]) * ndata
    
  }
  
  meanMat <- meanMat + cov(t(corsMat)) / nsims0
  
}

offdiags <- as.vector(meanMat[row(meanMat) != col(meanMat)])
diags <- diag(meanMat)

hist(offdiags * nsims0)
abline(v = offdVar(rhos), col = "red")
abline(v = mean(offdiags) * nsims0, col = "green")

hist(diags * nsims0)
abline(v = diagVar(rhos), col = "red")
abline(v = mean(diags) * nsims0, col = "green")

pvals <- pchisq(cr2sum, df = estvar, lower.tail = FALSE)
plot(-log10(seq_along(pvals) / (length(pvals) + 1)), -log10(sort(pvals)))
abline(0, 1, col = "red")
