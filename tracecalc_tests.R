library(rbenchmark)
source("mvrnormR.R")
source("tracecalcs.R")

m <- 500
mY <- 2
rho <- 0.5
nsims <- 1000
ndata <- 100
Beta <- 0
s2 <- 1

Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

set.seed(12345)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  # Data generation
  X <- mvrnormR(ndata, rep(0, m), Sig)
  Y <- Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
  
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
  
  Y4ColSums <- sum(Y4)
  X4ColSums <- colSums(X4)
  X1RowSums <- rowSums(X)
  X2RowSums <- rowSums(X2)
  
  # Trace calcs
  trace_uni()
  trace_uni_fast()
  
}

pvals <- pchisq(ndata * cors / as, df = bs, lower.tail = FALSE)
plot(-log10(seq_along(pvals) / (length(pvals) + 1)), -log10(sort(pvals)))
abline(0, 1)



corscov <- cov(t(corsmat)) * ndata
hist(diag(corscov))
abline(v = diagVar, col = "red")
hist(corscov[row(corscov) != col(corscov)])
abline(v = offdVar, col = "red")


Sig <- matrix(0, m, m)
diag(Sig) <- diagVar
Sig[row(Sig) != col(Sig)] <- offdVar
trSig <- sum(diag(Sig))
trSig2 <- sum(diag(crossprod(Sig)))
a <- trSig2 / trSig; b <- trSig^2 / trSig2

pvals <- pchisq(cors / a, df = b, lower.tail = FALSE)
#pvals <- pchisq(ndata * cors, df = m, lower.tail = FALSE)
hist(pvals)
#hist(sqrt(ndata) * colSums(corsmat) / sqrt(sum(Sig)))
