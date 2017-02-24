source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVars.R")
source("stdize.R")

m <- 20
rho <- 0
nsims <- 1000
ndata <- 1000
Beta <- 0
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- vars <- rep(0, nsims)
cormeans <- CM_vars <- rep(0, nsims)

set.seed(12345)

my_timers <- CM_timers <- rep(0, nsims)

for (sim in 1:nsims) {
  
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + rnorm(ndata, sd = sqrt(s2))
  
  
  # variance calculations
  corsums[sim] <- sum(cor(Yvars, Data))
  timer <- get_nanotime()
  vars[sim] <- varcalc1(Yvars, Data)
  my_timers[sim] <- get_nanotime() - timer
  cormeans[sim] <- mean(cor(Yvars, Data))
  
  timer <- get_nanotime()
  Ui <- (Yvars - mean(Yvars)) / sqrt(sum(Yvars^2))
  M_A <- stdize(t(Data))
  CM_vars[sim] <- makeVar(Ui, M_A)
  CM_timers[sim] <- get_nanotime() - timer
  
}

plot(CM_vars, vars)
