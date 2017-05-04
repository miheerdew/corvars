library(microbenchmark)
library(ggplot2)
m <- 100
my <- 100
rho <- 0.5
nsims <- 100
ndatas <- c(100)
Beta <- 1
s2 <- 1

source("mvrnormR.R")
source("varcalcs.R")
source("moment_calcs.R")
source("bh_reject.R")


SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
s2 <- s2 * 250 * (1 - rho + rho * 250)

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
    
    # Computing full update
    stats1 <- colSums(cor(X[ , 1:100], X)^2)
    vars1 <- apply(X, 2, function (Y) unlist(varmat_uni(Y, X[ , 1:100])))
    pvals <- pchisq(ndata * stats1 / vars1[1, ], df = vars1[2, ], lower.tail = FALSE)
    firstcors <- cor(X[ , 1], X)
    fishers <- sqrt(ndata - 3) * atanh(firstcors)
    fisherps <- pnorm(fishers, lower.tail = FALSE)
    B_new <- bh_reject(fisherps, 0.05)
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
      
      if (do_cluster_thres) {
        B_new2 <- cluster_thres(zs)
        if (length(B_new) > length(B_new2)) B_new <- B_new2
      }
      
      if (jaccard(B_old, B_new) == 0)
        break
      itCount <- itCount + 1
      
    }
    
    # Storing overlaps
    overlaps0[i, sim, 1] <- length(intersect(B_new, 1:100)) / length(B_new)
    overlaps0[i, sim, 2] <- length(intersect(B_new, 101:200)) / length(B_new)
    overlaps0[i, sim, 3] <- length(intersect(B_new, 201:300)) / length(B_new)
    
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
    
    # Completing first update
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
      B_new2 <- cluster_thres(zs)
      if (length(B_new) > length(B_new2)) B_new <- B_new2
      
      if (jaccard(B_old, B_new) == 0)
        break
      itCount <- itCount + 1
      
    }
    
    # Storing overlaps
    overlaps2[i, sim, 1] <- length(intersect(B_new, 1:100)) / length(B_new)
    overlaps2[i, sim, 2] <- length(intersect(B_new, 101:200)) / length(B_new)
    overlaps2[i, sim, 3] <- length(intersect(B_new, 201:300)) / length(B_new)
    
    if (overlaps2[i, sim, 1] < 0.75) break
    
    # Storing some marginal zs
    zs_comm[i, sim] <- zs1[101]
    zs_back[i, sim] <- zs1[201]
    
    if (runCM) {
    
      # Running CM
      CMout <- CM(t(X), max.groups = Inf, max.iter = 100, max.time = Inf)
      
      
      # Storing overlaps
      for (j in 1:min(length(CMout$DC_sets), 4)) {
        B <- CMout$DC_sets[[j]]
        overlapsCM[i, sim, 1, j] <- length(intersect(B, 1:100)) / length(B)
        overlapsCM[i, sim, 2, j] <- length(intersect(B, 101:200)) / length(B)
        overlapsCM[i, sim, 3, j] <- length(intersect(B, 201:300)) / length(B)
      }
      
      # Simulating null samples
      Xnull <- mvrnormR(ndata, rep(0, 3 * m), diag(3 * m))
      DCMout <- DCM(t(X), t(Xnull), max.groups = Inf, max.iter = 100, max.time = Inf)
      
      # Storing overlaps
      for (j in 1:min(length(DCMout$DC_sets), 4)) {
        B <- DCMout$DC_sets[[j]]
        overlapsDCM[i, sim, 1, j] <- length(intersect(B, 1:100)) / length(B)
        overlapsDCM[i, sim, 2, j] <- length(intersect(B, 101:200)) / length(B)
        overlapsDCM[i, sim, 3, j] <- length(intersect(B, 201:300)) / length(B)
      }
      
    }
    
  }
  
}


if (runCM) {
  percData <- data.frame("Method" = character(2 * length(ndatas)),
                         "Sample.Size" = numeric(2 * length(ndatas)),
                         "Count" = numeric(2 * length(ndatas)),
                         stringsAsFactors = FALSE)
  
  
  png("OLplotCM.png")
  par(mfrow = c(length(ndatas), 2))
  # Plotting stickiness
  for (i in seq_along(ndatas)) {
    
    
    plot(overlaps[i, , 2], apply(overlapsCM[i, , 1:2, 1], 1, max),
         ylab = "Max C set match to A1", xlab = "% A2 in U(A1)",
         main = paste0("CM behavior, n = ", ndatas[i]))
    
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
    labs(y = "# Times 2 comms were sticky",
         title = "CM/DCM sticky counts, by sample size")
  ggsave("sticky_counts.png", p)

}
