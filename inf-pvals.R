source("mvrnormR.R")
library(doParallel)

nsims <- 100
rho <- 0
ndata <- 600
alpha <- 0.05
MAXIT <- 10

mvrnormR <- function(n, mu, sigma) {
  ncols <- ncol(sigma)
  mu <- rep(mu, each = n) ## not obliged to use a matrix (recycling)
  mu + matrix(rnorm(n * ncols), ncol = ncols) %*% chol(sigma)
}

bh_reject <- function (pvals, alpha, conserv = TRUE) {

  m <- length(pvals)
  
  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals)
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals)
  }
  
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

sim <- function(ns) {
  foreach(n=ns, .combine = 'c') %:%
   foreach(sim=1:nsims,
          .combine=function(...) mean(c(...), na.rm=TRUE),
          .inorder = FALSE,
          .multicombine = TRUE) %dopar% {
     Sigma <- diag(n) #Completely independent case
     X <- mvrnormR(ndata, rep(0, n), Sigma)
     R0 <- crossprod(X) / ndata
     diag(R0) <- 0
     
     pvals <- function(A){
       vars <- rep(length(A), n)
       vars[A] <-  length(A)-1
       zs <- sqrt(ndata) * rowSums(R0[,A, drop=FALSE]) / sqrt(vars)
       pnorm(zs, lower.tail = FALSE)
     }
     
     max.pair <- as.vector(arrayInd(which.max(R0), .dim=dim(R0)))
     
      A_old <- c()
      A <- max.pair
      j <- 0
      while (!setequal(A, A_old) && length(A) > 1 && j < MAXIT){
        A_old <- A
        A <- bh_reject(pvals(A_old), alpha)
        j <- j + 1
      }
     #cat(sprintf("Size of fixed point %d\n", length(A)))
     ifelse(j >= MAXIT, NA, length(A) > 1)
  }
}

cl <- makePSOCKcluster(detectCores()-1)
setDefaultCluster(cl)
clusterExport(NULL, c('rho','ndata','mvrnormR','nsims', 'bh_reject', 'alpha', 'MAXIT'))
registerDoParallel(cl)

RNGkind("L'Ecuyer-CMRG")
set.seed(12347)


x <- c(100, 200, 300, 500, 800, 1000)
y <- sim(x)
plot(x, y, xlab="Number of variables", ylab="Probability of fixed point")
stopCluster(cl)
