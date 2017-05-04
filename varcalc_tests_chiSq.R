source("mvrnormR.R")
library(microbenchmark)
source("varcalcs_chiSq.R")
source("moment_calcs.R")
source("varcalcs.R")

m <- 20
mY <- 1
rho <- 0.5
nsims <- 1000
ndata <- 1000
Beta <- 0
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

cors <- vars <- as <- bs <- rep(0, nsims)
corsmat <- matrix(0, m, nsims)

#set.seed(12345)

my_timers <- CM_timers <- rep(0, nsims)

diagVar <- pyyjj(rho) + pxy(rho)^2 * (pyyyy(rho) + 2 * pyyjj(rho) + mu1111(rho)) / 4 - 
  pxy(rho) * (pyyyj(rho) + pyjjj(rho))
offdVar <- pyyjk(rho) + pxy(rho)^2 * (pyyyy(rho) + 2 * pyyjj(rho) + mu1122(rho)) / 4 - 
  pxy(rho) * (pyyyj(rho) + pyjjk(rho))

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
  vars[sim] <- varcalc1(Yvars[ , 1], Data)
  matvar <- varmat_uni(Yvars, Data)
  as[sim] <- matvar[[1]]
  bs[sim] <- matvar[[2]]
  allvars <- varcalc1_multi(Yvars, Data)
  corsi <- cor(Yvars, Data)
  cors[sim] <- tcrossprod(corsi)
  corsmat[ , sim] <- corsi
  
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
