source("mvrnormR.R")
source("tracecalcs.R")
library(bmdupdate)
library(rbenchmark)

# ndatas
ndatas <- seq(50, 500, 50)

# ms
ms <- seq(50, 500, 50)

# reps
reps <- 3

m <- 100
mY <- 200
rho <- 0.8
nsims <- 100
ndata <- 100
Beta <- 0
s2 <- 1


# Timer mats
timer_mat <- matrix(0, ncol = length(ms), nrow = length(ndatas))
rownames(timer_mat) <- ndatas
colnames(timer_mat) <- ms
timer_mats <- rep(list(timer_mat), 3)


for (i in seq_along(ndatas)) {
  
  for (j in seq_along(ms)) {
    
    ndata <- ndatas[i]
    m <- ms[j]
    
    cat("ndata =", ndata, "and m =", m, "\n")

    Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)
    
    # Data generation
    set.seed(12345)
    X <- mvrnormR(ndata, rep(0, m), Sig)
    Y <-
      Beta * rowSums(X) + matrix(rnorm(ndata * mY, sd = sqrt(s2)), ncol = mY)
    dx <- ncol(X)
    dy <- ncol(Y)
    
    X <- as.matrix(scale(X))
    Y <- as.matrix(scale(Y))
     
    # General Calcs
    n <- nrow(X)
    
    X2 <- X ^ 2
    X3 <- X ^ 3
    X4 <- X ^ 4
    Y2 <- Y ^ 2
    Y3 <- Y ^ 3
    Y4 <- Y ^ 4
    
    
    tX <- t(X)
    tY <- t(Y)
    tX2 <- t(X2)
    tX3 <- t(X3)
    XXt <- tcrossprod(X)
    XXt2 <- XXt ^ 2
    
    Y4ColSums <- colSums(Y4)
    X4ColSums <- colSums(X4)
    X2RowSums <- rowSums(X2)
    
    
    pvals0_fun <- function (A) {
      allr <- crossprod(Y, X[ , A, drop = FALSE])
      allrSums <- rowSums(allr)
      allr22 <- crossprod(Y2, X2[ , A, drop = FALSE])
      allr31 <- crossprod(Y3, X[ , A, drop = FALSE])
      source("tracecalcs.R", local = TRUE)
      trs1 <- sapply(1:dy, trace_indx)
      as <- trs1[2, ] / trs1[1, ]
      bs <- trs1[1, ] ^ 2 / trs1[2, ]
      cors <- crossprod(Y, X[ , A, drop = FALSE]) / (n - 1)
      pval <- pchisq(n * rowSums(cors^2) / as, df = bs, lower.tail = FALSE)
      return(pval)
    }
    
    bobj <- new(BmdUpdater, X, Y) 
    A <- 1:dx
    
    res <- benchmark(pvals0 = pvals0_fun(A),
                     pvals1 = bobj$pvals(A, FALSE),
                     pvals1p = bobj$pvals(A, TRUE),
                     replications = reps)
    
    timer_mats[[1]][i, j] <- res$elapsed[1]
    timer_mats[[2]][i, j] <- res$elapsed[2]
    timer_mats[[3]][i, j] <- res$elapsed[3]
    
  }
  
}

save(timer_mats, ndatas, ms, file = "bmdupdate_timers/bmdupdate_timer_data.RData")

# Plotting timers
library(reshape)
library(ggplot2)
library(grid)
library(gridExtra)

two2one <- melt(log2(timer_mats[[1]] / timer_mats[[2]]))
thr2two <- melt(log2(timer_mats[[2]] / timer_mats[[3]]))
names(two2one) <- names(thr2two) <- c("n", "m", "ratio")
maxplot <- max(abs(c(two2one$ratio, thr2two$ratio)))

p1 <- ggplot(two2one, aes(x = m, y = n, fill = ratio)) + geom_tile() + 
      geom_abline(slope = 1.2, intercept = 0, lwd = 3) + 
      guides(fill = guide_legend(title = "log2(meth1 / meth2)")) + 
      ggtitle("meth2 = BmdUpdater$pvals(par = F), meth1 = trace_indx") + 
      scale_fill_gradient(low = "blue", high = "red", limits = c(-maxplot, maxplot))
p2 <- ggplot(thr2two, aes(x = m, y = n, fill = ratio)) + geom_tile() + 
      geom_abline(slope = 1.2, intercept = 0, lwd = 3) + 
      guides(fill = guide_legend(title = "log2(meth1 / meth2)")) + 
      ggtitle("meth2 = BmdUpdater$pvals(par = T), meth1 = BmdUpdater$pvals(par = F)") + 
      scale_fill_gradient(low = "blue", high = "red", limits = c(-maxplot, maxplot))
plist <- list(p1, p2)

pdf("bmdupdate_timers/bmdupdate_timer_plot.pdf", width = 14)
grid.arrange(grobs = plist, ncol = 2, 
             top = textGrob("Timing Ratios", gp = gpar(fontsize = 20)))
dev.off()
