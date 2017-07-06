source("mvrnormR.R")
source("tracecalcs.R")

library(bmdupdate)
par <- TRUE

m <- 100
mY <- 1
rho <- 0.5
nsims <- 1000
ndata <- 100
Beta <- 4
s2 <- 1
#
Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
# 
pvals <- numeric(nsims)
set.seed(12345)
for (i in 1:nsims) {
  
  cat(i, "\n")
  
  # # Data generation
  X <- mvrnormR(ndata, rep(0, m), Sig)
  Y <-
    Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
  dx <- ncol(X)
  dy <- ncol(Y)
  
  bobj <- new(BmdUpdater, X, Y) 
  pvals[i] <- bobj$pvals(1:dx, par)
  
}

pquants <- c(1:nsims) / (nsims + 1)
plot(-log10(pquants), -log10(sort(pvals)))
abline(0, 1)
