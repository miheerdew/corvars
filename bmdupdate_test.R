source("mvrnormR.R")
source("tracecalcs.R")

m <- 100
mY <- 1
rho <- 0.9
nsims <- 100
ndata <- 100
Beta <- 0
s2 <- 1
# 
# 
Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
# 
# # Data generation
set.seed(12345)
X <- mvrnormR(ndata, rep(0, m), Sig)
Y <-
  Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
dx <- ncol(X)
dy <- ncol(Y)

# X <- as.matrix(scale(X))
# Y <- as.matrix(scale(Y))
# 
# # General Calcs
# n <- nrow(X)
# 
# X2 <- X ^ 2
# X3 <- X ^ 3
# X4 <- X ^ 4
# Y2 <- Y ^ 2
# Y3 <- Y ^ 3
# Y4 <- Y ^ 4
# 
# 
# tX <- t(X)
# tY <- t(Y)
# tX2 <- t(X2)
# tX3 <- t(X3)
# XXt <- tcrossprod(X)
# XXt2 <- XXt ^ 2
# 
# Y4ColSums <- colSums(Y4)
# X4ColSums <- colSums(X4)
# X2RowSums <- rowSums(X2)
# 
# allr <- crossprod(Y, X)
# allrSums <- rowSums(allr)
# allr22 <- crossprod(Y2, X2)
# allr31 <- crossprod(Y3, X)
# 
# 
# trs <- trace_indx(1)
# a <- trs[2] / trs[1]
# b <- trs[1] ^ 2 / trs[2]

library(bmdupdate)
par <- TRUE

bobj <- new(BmdUpdater, X, Y) 
A <- 1:dx #Test against the whole X set
bobj$pvals(A, par)
A <- dx + 1:dy #Test against the whole Y set
bobj$pvals(A, par)

# ??
