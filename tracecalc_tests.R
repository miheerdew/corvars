library(rbenchmark)
source("mvrnormR.R")
source("tracecalcs.R")

m <- 500
mY <- 1
rho <- 0.4
nsims <- 1000
ndata <- 100
Beta <- 0
s2 <- 1

Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

set.seed(12345)

as <- bs <- Ts <- ps <- numeric(nsims)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  # Data generation
  X <- mvrnormR(ndata, rep(0, m), Sig)
  Y <- Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
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
  allr22 <- crossprod(Y2, X2)
  allr31 <- crossprod(Y3, X)
  
  # Trace calcs
  #trace_uni(1)
  #mlist <- trace_uni_mlist(Y[ , 1], X)
  #unname(unlist(mlist[c("tr1", "tr2")]))
  #trace_large_x(Y[ , 1], X)
  #trace_kosher(Y[ , 1], X)
  
  trs <- trace_large_x_indx(1)
  as[sim] <- trs[2] / trs[1]
  bs[sim] <- trs[1]^2 / trs[2]
  Ts[sim] <- rowSums(allr^2 / (ndata - 1)^2)
  ps[sim] <- pchisq(ndata * Ts[sim] / as[sim], df = bs[sim], lower.tail = FALSE)
  
  if (sim == 1) {
    cat("doing timing for the first sim\n")
    timeres <- benchmark(trace_uni = trace_uni(1),
                         trace_uni_mlist = trace_uni_mlist(Y[ , 1], X),
                         trace_large_x = trace_large_x(Y[ , 1], X),
                         trace_large_x_indx = trace_large_x_indx(1))
  }
                     
  
}

plot(-log10(seq_along(ps) / (length(ps) + 1)), -log10(sort(ps)))
abline(0, 1)

