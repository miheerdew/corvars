source("mvrnormR.R")
source("tracecalcs.R")

set.seed(12345678)

n <- 500
nY <- 1000
rho <- 0.5
m <- 50
Beta <- 0
s2 <- 1
Sig <- diag(1 - rho, n) + matrix(rep(rho, n^2), ncol = n)

# Data generation

X <- mvrnormR(m, rep(0, n), Sig)
Y <- Beta * rowSums(X) + matrix(rnorm(m * nY, sd = sqrt(s2)), ncol = nY)

p <- profvis::profvis({
  m <- nrow(X)
  XXt <- tcrossprod(X)
  XXt2 <- XXt^2
  X2 <- X^2
  X3 <- X^3
  X4ColSums <- colSums(X^4)
  X2RowSums <- rowSums(X2)
  
  res <- sapply(1:nY, function(i) 
    trace_large_x(Y[,i], X, m, XXt, XXt2, X2, X3, X4ColSums, X2RowSums)
    )
})

#show(p)