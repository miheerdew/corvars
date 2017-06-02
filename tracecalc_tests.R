library(rbenchmark)
library(doParallel)
library(doRNG)
library(FNN)
source("mvrnormR.R")
source("tracecalcs.R")

#m <- 100
#mY <- 1
#rho <- 0.9
#nsims <- 100
#ndata <- 100
#Beta <- 0
#s2 <- 1


dist_to_uniform01 <- function(...){
  ks.test(c(...),"punif")$statistic
}

pvals <- function(m, rho, nsims, ndata, mY=1, s2=1, Beta=0){
  ps <- foreach(
    sim = 1:nsims,
    .combine = c,
    .export=c('mvrnormR', 'trace_indx', 'trace_uni_fast', 'trace_large_x_indx'),
    .packages = c('emulator')
  ) %dorng% {
    
    Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
    
    # Data generation
    X <- mvrnormR(ndata, rep(0, m), Sig)
    Y <-
      Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
    X <- as.matrix(scale(X))
    Y <- as.matrix(scale(Y))
    
    # General Calcs
    n <- nrow(X)
    
    X2 <- X ^ 2
    X3 <- X ^ 3
    X4 <- X ^ 4
    Y2 <- Y ^ 2
    Y3 <- Y ^ 3
    Y4 <- Y ^ 4
    
    
    tX <- t(X)
    tY <- t(Y)
    tX2 <- t(X2)
    tX3 <- t(X3)
    XXt <- tcrossprod(X)
    XXt2 <- XXt ^ 2
    
    Y4ColSums <- colSums(Y4)
    X4ColSums <- colSums(X4)
    X2RowSums <- rowSums(X2)
    
    allr <- crossprod(Y, X)
    allrSums <- rowSums(allr)
    allr22 <- crossprod(Y2, X2)
    allr31 <- crossprod(Y3, X)
    
    
    trs <- trace_indx(1)
    a <- trs[2] / trs[1]
    b <- trs[1] ^ 2 / trs[2]
    Tstat <- rowSums(allr ^ 2 / (ndata - 1) ^ 2)
    pchisq(ndata * Tstat / a, df = b, lower.tail = FALSE)
  }
  #attr(d, 'rng') <- NULL
  #d
  attr(ps, 'rng') <- NULL
  ps
}

setup <- function(cores=detectCores()-1){
  cl <<- makeCluster(cores, output="tracetest.txt")
  registerDoParallel(cl)
}

simple_test <- function(m, rho, ndata=1000, nsims=1000, plot=TRUE){
  set.seed(12345)
  ps <- pvals(m, rho, nsims, ndata)
  if(plot){
    hist(ps, freq=FALSE)
  }
  dist_to_uniform01(ps)
}



setup()

ndatas <- c(200, 500, 800, 1100)
rhos <- seq(0,0.9,by=0.1)
ms <-  c(1, 10, 50, 100, 200, 400)
 
N <- length(ms) * length(rhos) * length(ndatas)
simres <- data.frame(m=rep(NA,N), rho=rep(NA,N), ndata=rep(NA,N), dist=rep(NA,N), stringsAsFactors = FALSE)
i <- 1
 
for(m in ms){
  for(rho in rhos){
    for(ndata in ndatas){
      cat(sprintf("ndata=%d, rho=%f, m=%d\n", ndata, rho, m))
      simres[i,] <- list(m, rho, ndata, simple_test(m, rho, ndata))
      i <- i + 1
    }
  }
}

stopCluster(cl)

simres$ndata <- as.factor(simres$ndata)
simres$m <- as.factor(simres$m)
ggplot(simres, aes(y=dist, x=rho)) + geom_point(aes(color=ndata)) + facet_wrap( ~ m, ncol=3) + labs(title="Uniformity of ChiSqApprox pvalues", y="Distance from uniform", x="intra X correlations")


change_rhos <- function(rhos, m=100, ndata=1000, nsims=1000){
  res <- sapply(rhos, function(rho) simple_test(m, rho, ndata))
  plot(rhos, res)
  res
}
#Seems monotone in rho;