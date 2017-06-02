library(rbenchmark)
source("mvrnormR.R")
source("tracecalcs.R")
source("pointer_funs.R")
source("tracecalcs_pointers.R")

m <- 50
mY <- 25
rho <- 0.4
nsims <- 1000
ndata <- 10
Beta <- 0
s2 <- 1

Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

set.seed(12345)

as <- bs <- Ts <- ps <- numeric(nsims)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  # Data generation
  X0 <- mvrnormR(ndata, rep(0, m), Sig)
  Y0 <- Beta * rowSums(X0) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
  X0 <- newPointer(as.matrix(scale(X0)))
  Y0 <- newPointer(as.matrix(scale(Y0)))
  
  dx <- ncol(X0$value)
  dy <- ncol(Y0$value)
  
  # General Calcs
  n <- nrow(X0$value)
  
  X02 <- newPointer(X0$value^2)
  X03 <- newPointer(X0$value^3)
  X04 <- newPointer(X0$value^4)
  Y02 <- newPointer(Y0$value^2)
  Y03 <- newPointer(Y0$value^3)
  Y04 <- newPointer(Y0$value^4)
  
  tX0 <- newPointer(t(X0$value))
  tY0 <- newPointer(t(Y0$value))
  tX02 <- newPointer(t(X02$value))
  tY02 <- newPointer(t(Y02$value))
  tX03 <- newPointer(t(X03$value))
  tY03 <- newPointer(t(Y03$value))
  
  
  Y04ColSums <- newPointer(colSums(Y04$value))
  X04ColSums <- newPointer(colSums(X04$value))
  
  xyCors <- newPointer(crossprod(Y0$value, X0$value))
  txyCors <- newPointer(t(xyCors$value))
  xyCors22 <- newPointer(crossprod(Y02$value, X0$value))
  txyCors22 <- newPointer(crossprod(X02$value, Y0$value))
  xyCors31_y <- newPointer(crossprod(Y03$value, X0$value))
  xyCors31_x <- newPointer(crossprod(X03$value, Y0$value))
  
  # Choosing S
  S <- 1:50
  
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

