library(microbenchmark)
library(RCurl)

# Sourcing DCM functions
source("source_https.R")
repoDir <- "https://raw.githubusercontent.com/kbodwin/Differential-Correlation-Mining/master/DCM/R"
fNames <- c("CM", "prepData_CM", "init_CM", "run_CM", "resid_CM",
            "sanitize_CM", "stdize", "fisher", "makeVars", "emfa_DCM",
            "init_DCM", "resid_DCM", "run_DCM", "DCM", "sanitize_DCM",
            "prepData_DCM")
sNames <- file.path(repoDir, paste0(fNames, ".R"))
sourceOut <- source_https(sNames)

source("mvrnormR.R")
source("varcalcs.R")
source("moment_calcs.R")
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
overlapsCM <- array(0, dim = c(length(ndatas), nsims, 3, 4))
overlapsDCM <- array(0, dim = c(length(ndatas), nsims, 3, 4))

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
    
    # Running CM
    CMout <- CM(t(X), max.groups = Inf, max.iter = 100, max.time = Inf)
    
    
    # Storing overlaps
    for (j in 1:min(length(CMout$DC_sets), 4)) {
      B <- CMout$DC_sets[[j]]
      overlapsCM[i, sim, 1, j] <- length(intersect(B, 1:100)) / length(B_new)
      overlapsCM[i, sim, 2, j] <- length(intersect(B, 101:200)) / length(B_new)
      overlapsCM[i, sim, 3, j] <- length(intersect(B, 201:300)) / length(B_new)
    }
    
    # Simulating null samples
    Xnull <- mvrnormR(ndata, rep(0, 3 * m), diag(3 * m))
    DCMout <- DCM(t(X), t(Xnull), max.groups = Inf, max.iter = 100, max.time = Inf)
    
    # Storing overlaps
    for (j in 1:min(length(DCMout$DC_sets), 4)) {
      B <- DCMout$DC_sets[[j]]
      overlapsDCM[i, sim, 1, j] <- length(intersect(B, 1:100)) / length(B_new)
      overlapsDCM[i, sim, 2, j] <- length(intersect(B, 101:200)) / length(B_new)
      overlapsDCM[i, sim, 3, j] <- length(intersect(B, 201:300)) / length(B_new)
    }
    
    
    
  }
  
}


percData <- data.frame("Method" = character(6),
                       "Sample.Size" = numeric(6),
                       "Count" = numeric(6),
                       stringsAsFactors = FALSE)


png("OLplotCM.png")
par(mfrow = c(2, 3))
# Plotting stickiness
for (i in seq_along(ndatas)) {
  
  
  plot(overlaps[i, , 2], apply(overlapsCM[i, , 1:2, 1], 1, max),
       ylab = "Max C set match to A1", xlab = "% A2 in U(A1)",
       main = paste0("DCM behavior, n = ", ndatas[i]))
  
  percData[2 * (i - 1) + 1, "Method"] <- "CM"
  percData[2 * (i - 1) + 1, "Sample.Size"] <- ndatas[i]
  percData[2 * (i - 1) + 1, "Count"] <- 
    sum(apply(overlapsCM[i, , 1:2, 1], 1, max) <= .75)
  
  plot(overlaps[i, , 2], apply(overlapsDCM[i, , 1:2, 1], 1, max),
       ylab = "Max DC set match to A1", xlab = "% A2 in U(A1)",
       main = paste0("DCM behavior, n = ", ndatas[i]))
  
  percData[2 * (i - 1) + 2, "Method"] <- "DCM"
  percData[2 * (i - 1) + 2, "Sample.Size"] <- ndatas[i]
  percData[2 * (i - 1) + 2, "Count"] <- 
    sum(apply(overlapsDCM[i, , 1:2, 1], 1, max) <= .75)
    
}

dev.off()

p <- ggplot(percData, aes(x = Method, y = Count)) + 
  geom_bar(stat = "identity") + facet_grid(~Sample.Size) + 
  labs(y = "# Times 2 comms were sticky")
p

