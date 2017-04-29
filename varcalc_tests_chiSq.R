source("mvrnormR.R")
library(microbenchmark)
source("varcalcs_chiSq.R")
source("moment_calcs.R")

m <- 20
rho <- 0.5
nsims <- 1000
ndata <- 1000
Beta <- 0
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

cors <- vars <- rep(0, nsims)
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
  Yvars <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  vars[sim] <- varcalc1(Yvars, Data)
  corsi <- cor(Yvars, Data)
  cors[sim] <- tcrossprod(corsi)
  corsmat[ , sim] <- corsi
  
}

corscov <- cov(t(corsmat)) * ndata
hist(diag(corscov))
abline(v = diagVar, col = "red")
hist(corscov[row(corscov) != col(corscov)])
abline(v = offdVar, col = "red")



#-------------------------------------------------------------------------------------------------------

cors0 <- cors
vars0 <- vars


m <- 20
rho <- 0.5
nsims <- 1000
ndata <- 4000
Beta <- 0
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

cors <- vars <- rep(0, nsims)

#set.seed(12345)

my_timers <- CM_timers <- rep(0, nsims)

pvar <- m * (pyyjj(rho) + pxy(rho)^2 * (pyyyy(rho) + 2 * pyyjj(rho) + pyyyy(rho)) / 4 - 
               pxy(rho) * (pyyyj(rho) + pyyyj(rho)))

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  vars[sim] <- varcalc1(Yvars, Data)
  cors[sim] <- tcrossprod(cor(Yvars, Data))
  
  
}

cors1 <- cors
vars1 <- vars

plot(sort(cors0), sort(cors1))
abline(0, 1)
plot(sort(vars0), sort(vars1))
abline(0, 1)
