source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVars.R")
source("stdize.R")
source("bh_reject.R")

m <- 100
my <- 100
rho <- 0.5
nsims <- 100
ndatas <- seq(1000, 3000, 1000)
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
s2 <- s2 * 250 * (1 - rho + rho * 250)

zs_comm <- zs_back <- matrix(0, ndatas, nsims)
overlaps <- array(0, dim = c(length(ndatas), nsims, 3))

set.seed(12345)

for (i in seq_along(ndatas)) {
  
  ndata <- ndatas[i]
  
  cat("ndata =", ndata, "\n")

  for (sim in 1:nsims) {
    
    if (sim %% 10 == 0)
      cat("----sim ", sim, "\n")
    
    Data1 <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
    #Yvars1 <- Beta * rowSums(Data1) + matrix(rnorm(ndata * my, sd = sqrt(s2)),
                                            #ncol = my)
    
    Data2 <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
    #Yvars2 <- Beta * rowSums(Data2) + matrix(rnorm(ndata * my, sd = sqrt(s2)),
                                            #ncol = my)
    
    Data3 <- mvrnormR(ndata, rep(0, m), diag(1, m))
    #Yvars3 <- matrix(rnorm(ndata * my, sd = sqrt(s2)), ncol = my)
    
    X <- cbind(Data1, Data2, Data3)
    #Y <- cbind(Yvars1, Yvars2, Yvars3)
    
    # Computing update
    Us <- X
    MA <- X[ , 1:100]
    corsums1 <- rowSums(cor(Us, MA))
    vars1 <- varcalc1_multi(Us, MA)
    zs1 <- sqrt(ndata) * (corsums1) / sqrt(vars1)
    ps1 <- pnorm(zs1, lower.tail = FALSE)
    B1 <- bh_reject(ps1, 0.05)
    
    # Storing overlaps
    overlaps[i, sim, 1] <- length(intersect(B1, 1:100)) / length(B1)
    overlaps[i, sim, 2] <- length(intersect(B1, 101:200)) / length(B1)
    overlaps[i, sim, 3] <- length(intersect(B1, 201:300)) / length(B1)
    
    # Storing some marginal zs
    zs_comm[i, sim] <- zs1[101]
    zs_back[i, sim] <- zs1[201]
    
  }
  
}

