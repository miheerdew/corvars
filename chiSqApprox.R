source("mvrnormR.R")

nsims <- 1000
m <- 10
rho <- 0.4
Sig <- diag(1 - rho, m) + matrix(rep(rho, m^2), ncol = m)

Data <- mvrnormR(nsims, mu = rep(0, m), sigma = Sig)
stats <- rowSums(Data^2)
Varcalc <- m * (3 + (m - 1) * (1 + 2 * rho^2)) - m^2
var(stats)
Varcalc

trSig <- sum(diag(Sig))
trSig2 <- sum(diag(crossprod(Sig)))
a <- trSig2 / trSig; b <- trSig^2 / trSig2

#compx <- rchisq(nsims, df = Varcalc / 2)
#compx <- rgamma(nsims, shape = m^2 / Varcalc, rate = m / Varcalc)
compx <- a * rchisq(nsims, df = b)
mean(compx)
mean(stats)

if (TRUE) {
  #pvals <- pgamma(stats, shape = m^2 / Varcalc, rate = m / Varcalc, lower.tail = FALSE)
  #pvals <- pexp(stats, rate = 1 / sqrt(Varcalc), lower.tail = FALSE)
  pvals <- pchisq(stats / a, df = b, lower.tail = FALSE)
} else {
  pvals <- pnorm(stats, mean = m, sd = sqrt(Varcalc), lower.tail = FALSE)
}
plot(-log10(seq_along(pvals) / (length(pvals) + 1)), -log10(sort(pvals)))
abline(0, 1)
