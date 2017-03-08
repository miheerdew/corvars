source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVars.R")
source("stdize.R")

m <- 20
my <- 10
rho <- 0
nsims <- 1000
ndata <- 1000
Beta <- 0
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

corsums <- vars <- rep(0, nsims * my)
cormeans <- CM_vars <- rep(0, nsims * my)

set.seed(12345)

my_timers <- CM_timers <- rep(0, nsims)

for (sim in 1:nsims) {
  
  posindx <- (sim - 1) * my + 1:my
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + matrix(rnorm(ndata * my, sd = sqrt(s2)),
                                         ncol = my)
  
  
  # variance calculations
  corsums[posindx] <- rowSums(cor(Yvars, Data))
  timer <- get_nanotime()
  vars[posindx] <- varcalc1_multi(Yvars, Data)
  my_timers[sim] <- get_nanotime() - timer

  
  cormeans[posindx] <- rowMeans(cor(Yvars, Data))  
  timer <- get_nanotime()
  Ui <- stdize(t(Yvars))
  M_A <- stdize(t(Data))
  CM_vars[posindx] <- makeVars(Ui, M_A)
  CM_timers[sim] <- get_nanotime() - timer
  
}

zs <- sqrt(ndata) * (corsums - m * pxy(rho)) / sqrt(vars)
CM_zs <- cormeans / sqrt(CM_vars)

plot(CM_zs, zs)

sum((CM_zs - zs)^2)

