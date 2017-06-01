source("mvrnormR.R")
source("tracecalcs.R")
library(rbenchmark)

set.seed(12345678)

#Breaking point ndata > 1.2*m
m <- 200
nY <- 1
rho <- 0.5
ndata <- 200
Beta <- 0.2
s2 <- 1
Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

# Data generation

X <- mvrnormR(ndata, rep(0, m), Sig)
Y <- Beta * rowSums(X) + matrix(rnorm(ndata * nY, sd = sqrt(s2)), ncol = nY)


  XXt <- tcrossprod(X)
  XXt2 <- XXt^2
  X2 <- X^2
  X3 <- X^3
  X4ColSums <- colSums(X^4)
  X2RowSums <- rowSums(X2)
  tX <- t(X)

  X <- mvrnormR(ndata, rep(0, m), Sig)
  Y <- Beta * rowSums(X) + matrix(rnorm(ndata * nY, sd = sqrt(s2)), ncol = nY)
  X <- as.matrix(scale(X))
  Y <- as.matrix(scale(Y))
  
  # General Calcs
  n <- nrow(X)
  
  X2 <- X^2
  X3 <- X^3
  X4 <- X^4
  Y2 <- Y^2
  Y3 <- Y^3
  Y4 <- Y^4
  
  tX <- t(X)
  tY <- t(Y)
  tX2 <- t(X2)
  tX3 <- t(X3)
  XXt <- tcrossprod(X)
  XXt2 <- XXt^2
  
  Y4ColSums <- colSums(Y4)
  X4ColSums <- colSums(X4)
  X1RowSums <- rowSums(X)
  X2RowSums <- rowSums(X2)
  
  allr <- crossprod(Y, X)
  allrSums <- rowSums(allr)
  allr22 <- crossprod(Y2, X2)
  allr31 <- crossprod(Y3, X)
  
#  res <- sapply(1:nY, function(i)
#    trace_large_x(Y[,i], X, m, XXt, XXt2, X2, X3, X4ColSums, X2RowSums, tX)
#    )
#p <- profvis::profvis({
#  trace_uni_fast(1)
#})


assertthat::are_equal(trace_uni(1), trace_uni_fast(1))
print(benchmark(large_x=trace_large_x_indx(1), 
                uni_fast=trace_uni_fast(1), 
#                uni=trace_uni(1),
                replications=10))
#print(p)
