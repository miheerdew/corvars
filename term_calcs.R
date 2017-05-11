source("mvrnormR.R")

# Parameters
n <- 500
nY <- 2
rho <- 0.5
m <- 100
Beta <- 0
s2 <- 1
Sig <- diag(1 - rho, n) + matrix(rep(rho, n^2), ncol = n)

# Data generation
X <- mvrnormR(ndata, rep(0, n), Sig)
Y <- Beta * rowSums(X) + matrix(rnorm(m * nY, sd = sqrt(s2)), ncol = nY)

# General Calcs
X <- scale(X)
Y <- scale(Y)

X2 <- X^2
X3 <- X^3
X4 <- X^4
Y2 <- Y^2
Y3 <- Y^3
Y4 <- Y^4


tX <- t(X)
tY <- t(Y)
tX2 <- t(X2)
tX3 <- t(X3)
XXt <- tcrossprod(X)
XXt2 <- XXt^2

Y4ColSums <- sum(Y4)
X4ColSums <- colSums(X4)
X1RowSums <- rowSums(X)
X2RowSums <- rowSums(X2)

# Cor calcs
xyCors <- crossprod(X, Y) / (m - 1)
xyCors2 <- xyCors^2
rX <- crossprod(tX, xyCors)
r2X <- crossprod(tX, xyCors2)
r2X2 <- crossprod(tX2, xyCors2)
rX3 <- crossprod(tX3, xyCors)


# The following is a sandbox to calculate various cross-terms and check them
source("tracecalcs.R")
rm(trace_uni_fast, trace_uni)

# This is the term we need to square (see the Google doc):

# Sjk = ruujk + ruj * ruk(ruuuu + ruujj + ruukk + rjjkk) / 4 -
#       ruj(ruuuk + rujjk) / 2 - ruk(ruuuj + rukkj) / 2

# I have written a function to calculate 9 |B| x |B| matrices,
# containing the components of the above sum, returned in a list.

# The list contains the elements A ... C4, which correspond to Sjk:
# S = A + (B1 + B2 + B3 + B4) / 4 - (C1 + C2) / 2 - (C3 + C4) / 2
# Recall that sum(S) = Var(R(Y, X)), where here Y is a single Y vector

# The usage of the function is:
mlist <- trace_uni_mlist(Y[ , 1], X)
# ...it does not take the full Y matrix.

# This function is useful to check the accuracy of cross-term calculation.
# As a first example, let's try to compute sum(A * A)

AA <- crossprod(Y2[ , 1], crossprod(XXt2, Y2[ , 1]))
sum(mlist$A^2)
AA

# The above is m^2 * n, as desired. Unfortunately, not all cross terms are as nice. 
# The cross-terms can be divided into 3 categories:
#   (1) Symmetric: AA = ruujk * ruujk
#   (2) Non-symmetric, half-pairwise: AB2 = ruujk * ruj * ruk * ruujj
#   (3) Non-symmetric, full-pairwise: AC2 = ruujk * ruj * rujjk

# Cross-terms of category (2) can be obtained by pre-computing the non-pairwise term.
# Cross-terms of category (3) are more difficult to deal with.
# I will attempt AC2. See the Google doc for the derivation.
AC2mat <- crossprod(xyCors[ , 1] * tX, tX2)
AC2 <- crossprod(Y[ , 1], crossprod(AB2mat * XXt, Y2[ , 1]))
AC2
sum(mlist$A * mlist$C2)

# The issue is that the calculation of AC2 mat is not easily vectorizable across
# multiple Y vectors. One solution would be to do an sapply over Y column indices, but
# that is probably not desirable.