source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVars.R")
source("stdize.R")

n <- 1000
m <- 10
nsims <- 100
rho <- 0
Sigma <- diag(1 - rho, m) + matrix(rho, m, m)


my_vars <- CM_vars <- teststats <- teststats2 <- rep(0, nsims)


for (i in 1:nsims) {
  
  Data <- mvrnormR(n, rep(0, m), Sigma)
  teststats[i] <- mean(cor(Data[ , 10], Data[ , 1:9]))
  
  my_vars[i] <- varcalc1(Data[ , 10], Data[ , 1:9])
  teststats2[i] <- sum(cor(Data[ , 10], Data[ , 1:9]))
  
  M_A <- stdize(t(Data[ , 1:9]))
  Ui <- Data[ , 10] / sqrt(sum(Data[ , 10]^2))
  
  CM_vars[i] <- makeVar(Ui, M_A)
  
}