library(rbenchmark)
source("mvrnormR.R")
source("tracecalcs.R")
source("pointer_funs.R")

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
  X <- newPointer(as.matrix(scale(X)))
  Y <- newPointer(as.matrix(scale(Y)))
  
  # General Calcs
  n <- nrow(X)
  
  X2 <- newPointer(X^2)
  X3 <- newPointer(X^3)
  X4 <- newPointer(X^4)
  Y2 <- newPointer(Y^2)
  Y3 <- newPointer(Y^3)
  Y4 <- newPointer(Y^4)
  
  
  tX <- newPointer(t(X))
  tY <- newPointer(t(Y))
  tX2 <- newPointer(t(X2)); tY2 <- newPointer(t(Y2))
  tX3 <- newPointer(t(X3)); tY3 <- newPointer(t(Y3))
  XXt <- newPointer(tcrossprod(X)); YYt <- newPointer(tcrossprod(Y))
  XXt2 <- newPointer(XXt^2); YYt2 <- newPointer(YYt^2)
  
  Y4ColSums <- newPointer(colSums(Y4))
  X4ColSums <- newPointer(colSums(X4))
  X1RowSums <- newPointer(rowSums(X)); X2RowSums <- newPointer(rowSums(X2))
  Y1RowSums <- newPointer(rowSums(Y)); Y2RowSums <- newPointer(rowSums(Y2))
  
  allr <- crossprod(Y, X)
  allr22 <- crossprod(Y2, X2)
  allr31 <- crossprod(Y3, X)
  allr13 <- crossprod(Y, X3)
  
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

