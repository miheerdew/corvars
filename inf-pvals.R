source("mvrnormR.R")
source("varcalcs.R")
library(parallel)

nsims <- 100
rho <- 0.1
ndata <- 600


pval <- function(y, X){
  ndata <- nrow(X)
  pnorm(sqrt(ndata) * sum(cor(y,X))/sqrt(varcalc1(y, X)), lower.tail = FALSE)
}

sim <- function(ns) {
  foreach(n=ns, .combine = 'c') %:%
   foreach(sim=1:nsims,
          .combine=function(...) mean(c(...)),
          .inorder = FALSE,
          .multicombine = TRUE) %dopar% {
     Sigma <- diag(1 - rho, n) + matrix(rep(rho, n^2), ncol = n)
     X <- mvrnormR(ndata, rep(0, n), Sigma)
     y <- rnorm(ndata)
     r <- as.vector(crossprod(y, X)/ndata)
     ids <- order(r, decreasing = TRUE)
    pval(y, X[,ids[1], drop=FALSE])
  }
}

cl <- makePSOCKcluster(detectCores()-1)
setDefaultCluster(cl)
clusterExport(NULL, c('rho','ndata','mvrnormR','nsims', 'pval','varcalc1'))
registerDoParallel(cl)

RNGkind("L'Ecuyer-CMRG")
set.seed(12347)

x <- c(seq(100, 900, by=100), seq(1000, 4000, by=500))
y <- sim(x)
plot(x, 1/y, xlab="n", ylab="1/min(p-value)")
stopCluster(cl)