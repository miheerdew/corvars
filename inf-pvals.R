source("mvrnormR.R")
source("varcalcs.R")
library(parallel)

nsims <- 60
rho <- 0.4
ndata <- 600


pval <- function(y, X){
  ndata <- nrow(X)
  pnorm(sqrt(ndata) * sum(cor(y,X))/sqrt(varcalc1(y, X)), lower.tail = FALSE)
}

sim <- function(n) {
  ps <- numeric(nsims)
  for(i in 1:nsims){
     Sigma <- diag(1 - rho, n) + matrix(rep(rho, n^2), ncol = n)
     X <- mvrnormR(ndata, rep(0, n), Sigma)
  #   y <- rnorm(ndata)
  #   r <- as.vector(crossprod(y, X)/ndata)
  #   ids <- order(r, decreasing = TRUE)
  #   ps[i] <- pval(y, X[,ids[1], drop=FALSE])
  }
  mean(ps)
  ndata
}

RNGkind("L'Ecuyer-CMRG")
set.seed(12346)

x <- seq(100, 200, by=100)
y <- unlist(mclapply(x, sim, mc.silent = FALSE))
plot(x, -log10(y), xlab="n", ylab="-log10(min-pval)")
print(y)