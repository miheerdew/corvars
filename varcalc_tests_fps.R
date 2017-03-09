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
JP <- FALSE

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
s2 <- s2 * 250 * (1 - rho + rho * 250)

zs_comm <- zs_back <- matrix(0, ndatas, nsims)
overlaps <- overlaps2 <- array(0, dim = c(length(ndatas), nsims, 3))

# Seed generations
set.seed(12345)
seeds <- lapply(seq_along(ndatas), 
                function (i) sample(1:1e7, nsims, replace = FALSE))

for (i in seq_along(ndatas)) {
  
  ndata <- ndatas[i]
  
  cat("ndata =", ndata, "\n")

  for (sim in 1:nsims) {
    
    if (sim %% 10 == 0)
      cat("----sim ", sim, "\n")
    
    set.seed(seeds[[i]][sim])
    
    Data1 <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
    
    Data2 <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
    
    Data3 <- mvrnormR(ndata, rep(0, m), diag(1, m))
    
    X <- cbind(Data1, Data2, Data3)
    
    # Computing first update
    Us <- X
    MA <- X[ , 1:100]
    if (JP) {
      corsums1 <- rowSums(cor(Us, MA))
      vars1 <- varcalc1_multi(Us, MA) / ndata
    } else {
      corsums1 <- rowMeans(cor(Us, MA))  
      Us <- stdize(t(Us))
      MA <- stdize(t(MA))
      vars1 <- makeVars(Us, MA)
    }
    zs1 <- corsums1 / sqrt(vars1)
    ps1 <- pnorm(zs1, lower.tail = FALSE)
    B1 <- bh_reject(ps1, 0.05)
    
    # Storing overlaps
    overlaps[i, sim, 1] <- length(intersect(B1, 1:100)) / length(B1)
    overlaps[i, sim, 2] <- length(intersect(B1, 101:200)) / length(B1)
    overlaps[i, sim, 3] <- length(intersect(B1, 201:300)) / length(B1)
    
    # Completing updates
    B_new <- B1
    itCount <- 1
    
    repeat {
      
      B_old <- B_new
      cat("--------iteration", itCount, "\n")

      # Computing update
      Us <- X
      MA <- X[ , B_new]
      if (JP) {
        corsums <- rowSums(cor(Us, MA))
        vars <- varcalc1_multi(Us, MA) / ndata
      } else {
        corsums <- rowMeans(cor(Us, MA))  
        Us <- stdize(t(Us))
        MA <- stdize(t(MA))
        vars <- makeVars(Us, MA)
      }
      zs <- corsums / sqrt(vars)
      ps <- pnorm(zs, lower.tail = FALSE)
      B_new <- bh_reject(ps, 0.05)
      
      if (jaccard(B_old, B_new) == 0)
        break
      itCount <- itCount + 1
      
    }
    
    # Storing overlaps
    overlaps2[i, sim, 1] <- length(intersect(B_new, 1:100)) / length(B_new)
    overlaps2[i, sim, 2] <- length(intersect(B_new, 101:200)) / length(B_new)
    overlaps2[i, sim, 3] <- length(intersect(B_new, 201:300)) / length(B_new)
    
    # Storing some marginal zs
    zs_comm[i, sim] <- zs1[101]
    zs_back[i, sim] <- zs1[201]
    
  }
  
}

