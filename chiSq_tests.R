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

Data1 <- mvrnormR(ndata, rep(0, bX), SigmaX)
cr2sum1 <- rowSums(Data1^2)
Data2 <- mvrnormR(ndata, rep(0, bX), diag(bX))
#cr2sum2 <- diag(Data2 %*% SigmaX %*% t(Data2))
#cr2sum2 <- diag(Data2 %*% diag(bX) %*% t(Data2))
cr2sum2 <- rowSums(Data2^2)

plot(sort(cr2sum1), sort(cr2sum2))
abline(0, 1)
