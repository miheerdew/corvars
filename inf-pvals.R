source("mvrnormR.R")
source("bh_reject.R")
library(doParallel)

nsims <- 50
rho <- 0
ndata <- 600
alpha <- 0.05

sim <- function(ns) {
  foreach(n=ns, .combine = 'c') %:%
   foreach(sim=1:nsims,
          .combine=function(...) mean(c(...)),
          .inorder = FALSE,
          .multicombine = TRUE) %dopar% {
     Sigma <- diag(n) #Completely independent case
     X <- mvrnormR(ndata, rep(0, n), Sigma)
     R0 <- crossprod(X) / ndata
     diag(R0) <- 0
     
     pvals <- function(A){
       vars <- rep(length(A), n)
       vars[A] <-  length(A)-1
       zs <- sqrt(ndata) * rowSums(R0[,A]) / sqrt(vars)
       pnorm(zs, lower.tail = FALSE)
     }
     
     max.pair <- as.vector(arrayInd(which.max(R0), .dim=dim(R0)))
     
     A_old <- c()
     A <- max.pair
     #j <- 0
     while (!setequal(A, A_old) && length(A) > 0){
       A_old <- A
       A <- bh_reject(pvals(A_old), alpha)
     }
     #cat(sprintf("Size of fixed point %d\n", length(A)))
     length(A) > 0
  }
}

cl <- makePSOCKcluster(detectCores()-1)
setDefaultCluster(cl)
clusterExport(NULL, c('rho','ndata','mvrnormR','nsims', 'bh_reject', 'alpha'))
registerDoParallel(cl)

RNGkind("L'Ecuyer-CMRG")
set.seed(12347)


x <- seq(100, 500, 100)#c(seq(100, 1000, by=100), seq(1000, 4000, by=500))
y <- sim(x)
plot(x, y, xlab="n", ylab="Probability of fixed point")
stopCluster(cl)