source("mvrnormR.R")
library(microbenchmark)
source("varcalcs.R")
source("moment_calcs.R")
source("makeVars.R")
source("stdize.R")

m <- 50
my <- 50
rho <- 0.5
nsims <- 1000
ndata <- 1000
Beta <- 1
s2 <- 1

SigmaX_hom <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

kbTimers <- matrix(0, nsims, 8)
jpTimers <- matrix(0, nsims, 8)

set.seed(12345)

for (sim in 1:nsims) {
  
  # Make the data
  cat("sim ", sim, "\n")
  
  Data <- mvrnormR(ndata, rep(0, m), SigmaX_hom)
  Yvars <- Beta * rowSums(Data) + matrix(rnorm(ndata * my, sd = sqrt(s2)),
                                         ncol = my)
  
  # KB calcs
  
  # Initial calcs
  timer <- get_nanotime()
  Uis <- stdize(t(Yvars))
  M_A <- stdize(t(Data))
  kbTimers[sim, 1] <- get_nanotime() - timer
  
  # Frontmatter
  timer <- get_nanotime()
  n1 = ncol(Uis)
  r1s = M_A %*% t(Uis)
  k = nrow(r1s)
  Uis = Uis*sqrt(n1)
  M_A = M_A*sqrt(n1)
  W = t(array(colMeans(M_A), c(ncol(Uis), nrow(Uis))))
  Y = (t(r1s) %*% M_A^2)/k
  rA = colMeans(r1s)
  kbTimers[sim, 2] <- get_nanotime() - timer

  # Calcs
  timer <- get_nanotime()
  a <- 1/4*rA^2*rowMeans(Uis^4)
  kbTimers[sim, 4] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  b <- rA*rowMeans(Y/2*Uis^2) 
  kbTimers[sim, 5] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  c <- rowMeans(W^2*Uis^2 + Y^2/4 - W*Y*Uis)
  kbTimers[sim, c(3, 6, 8)] <- (get_nanotime() - timer) / 3
  
  timer <- get_nanotime()
  d <- -rA*rowMeans(W*Uis^3)
  kbTimers[sim, 7] <- get_nanotime() - timer
  
  # JP calcs
  
  # Initial calcs
  timer <- get_nanotime()
  n <- nrow(Data)
  Xmat <- scale(Data)
  Ymat <- scale(Yvars)
  jpTimers[sim, 1] <- get_nanotime() - timer
  
  # Frontmatter
  timer <- get_nanotime()
  # General calcs
  xyCors <- cor(Ymat, Xmat) #r1s
  SxyCors <- rowSums(xyCors) #rA
  xRowSum <- rowSums(Xmat) #W
  xRowSum2 <- tcrossprod(xyCors, Xmat^2) #Y
  jpTimers[sim, 2] <- get_nanotime() - timer

  # Calcs
  timer <- get_nanotime()
  star1 <- crossprod(Ymat^2, xRowSum^2)
  jpTimers[sim, 3] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  star2 <- colSums(Ymat^4) * rowSums(xyCors)^2
  jpTimers[sim, 4] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  star3 <- 2 * xyCors * colSums(Ymat^2 * t(xRowSum2))
  jpTimers[sim, 5] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  star4 <- rowSums(xRowSum2^2)
  jpTimers[sim, 6] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  dagger1 <- rowSums(xyCors) * crossprod(Ymat^3, xRowSum)
  jpTimers[sim, 7] <- get_nanotime() - timer
  
  timer <- get_nanotime()
  dagger2 <- colSums(xRowSum * t(xRowSum2) * Ymat)
  jpTimers[sim, 8] <- get_nanotime() - timer
  
  
}

apply(kbTimers, 2, median) / 1e6

apply(jpTimers, 2, median) / 1e6
