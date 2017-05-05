trace_uni <- function (i) {
  
  # General calcs
  xyCors <- as.vector(crossprod(Y[ , i], X)) / (n - 1)
  r <- sum(xyCors)
  rX2 <- crossprod(tX2, xyCors)
  
  # star 1
  S1 <- crossprod(as.vector(Y[ , i]) * X)
  
  # star 2
  S2 <- Y4ColSums * tcrossprod(xyCors)
  
  # star 3
  S3 <- tcrossprod(xyCors, xyCors * colSums(X2 * as.vector(Y2[ , i])))
  
  # star 4
  S4 <- tcrossprod(xyCors * t(X2))
  
  # dagger 1
  D1 <- tcrossprod(xyCors, colSums(X * as.vector(Y3[ , i])))
  
  # dagger 2
  D2 <- t(crossprod(as.vector(Y[ , i]) * X, X2)) * xyCors
  
  bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)
  
  trSig <- sum(diag(bigmat))
  trSig2 <- sum(diag(crossprod(bigmat)))
  
  return(c(trSig, trSig2))
  
}

trace_uni_mlist <- function (Yvec, Xmat) {
  
  Xmat <- as.matrix(scale(Xmat))
  Yvec <- as.vector(scale(Yvec))
  
  n <- nrow(Xmat)
  
  if (n != length(Yvec))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- as.vector(cor(Yvec, Xmat))
  y4 <- sum(Yvec^4)
  xRowSum <- rowSums(Xmat)
  xRowSum2 <- tcrossprod(xyCors, Xmat^2)
  
  # star 1
  S1 <- crossprod(Yvec * Xmat)
  
  # star 2
  S2 <- y4 * tcrossprod(xyCors)
  
  # star 3
  S3 <- tcrossprod(xyCors, xyCors * colSums(Xmat^2 * Yvec^2))
  
  # star 4
  S4 <- tcrossprod(xyCors * t(Xmat^2))
  
  # dagger 1
  D1 <- tcrossprod(xyCors, colSums(Xmat * Yvec^3))
  
  # dagger 2
  D2 <- t(crossprod(Yvec * Xmat, Xmat^2)) * xyCors
  
  bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)
  
  return(list(A = S1, 
              B1 = S2 / 4, B2 = S3 / 4, B3 = t(S3) / 4, B4 = S4 / 4,
              C1 = -t(D1) / 2, C2 = -D1 / 2, C3 = -t(D2) / 2, C4 = -D2 / 2))
  
}

trace_uni_fast <- function () {
  
  # General calcs
  xyCors <- crossprod(X, Y) / (n - 1)
  xyCors2 <- xyCors^2
  rX <- crossprod(tX, xyCors)
  r2X <- crossprod(tX, xyCors2)
  r2X2 <- crossprod(tX2, xyCors2)
  rX3 <- crossprod(tX3, xyCors)
  
  
  #AA
  AA <- crossprod(Y2, crossprod(XXt2, Y2))
  
  
  
  #AB1
  AB1 <- Y4ColSums * crossprod(Y2, rX^2) / 4
  
  
  
  #AB2
  AB2array <- do.call(cbind, lapply(1:n, function (i) tX[ , i] * tX2))
  #AB2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
  AB2mat <- tcrossprod(
  AB2list <- lapply(1:ncol(AB2mat), function (i) matrix(AB2mat[ , i], nrow = n, byrow = TRUE))
  AB2 <- crossprod(Y2, crossprod(AB2mat, Y2 * rX)) / 4
  
  #AB3
  AB3 <- AB2
  
  
  
  #AB4
  AB4 <- crossprod(rep(1, n), crossprod(AB2mat^2, Y2)) / 4
  
  
  
  #AC1
  AC1vec1 <- crossprod(Y2 * rX, X)
  AC1vec2 <- crossprod(X, Y3)
  AC1 <- -AC1vec1 %*% AC1vec2 / 2
  
  #AC2
  AC2 <- AC1
  
  
  
  #AC3
  AC3 <- -crossprod(Y, crossprod(AB2mat * XXt, Y2)) / 2
  
  
  
  #AC4
  AC4 <- AC3
  
  
  
  #B1^2
  B1B1 <- Y4ColSums^2 * sum(xyCors2)^2 / 16
  
  
  
  #B1B2
  B1B2 <- Y4ColSums * crossprod(Y2, r2X2) * sum(xyCors2) / 16
  
  
  
  #B1B3
  B1B3 <- B1B2
  
  
  
  #B1B4
  B1B4 <- Y4ColSums * crossprod(r2X2) / 16
  
  
  
  #B1C1
  B1C1 <- -Y4ColSums * sum(xyCors2) * crossprod(Y3, rX) / 8
  
  

  #B1C2
  B1C2 <- B1C1
  
  
  
  #B1C3
  B1C3 <- -Y4ColSums * crossprod(Y, rX * r2X2) / 8
  
  
  
  #B1C4
  B1C4 <- B1C3
  
  
  
  #B2^2
  B2B2array <- do.call(cbind, lapply(1:n, function (i) tX2[ , i] * tX2))
  B2B2mat <- matrix(crossprod(B2B2array, xyCors2), n, byrow = TRUE)
  B2B2 <- sum(xyCors2) * crossprod(Y2, crossprod(B2B2mat, Y2)) / 16
  
  #B2B3
  B2B3 <- crossprod(Y2, r2X2)^2 / 16
  
  
  
  #B2B4
  B2B4 <- crossprod(Y2, crossprod(B2B2mat, r2X2)) / 16
  
  
  
  #B2C1
  B2C1 <- -crossprod(Y3, rX) * crossprod(Y2, r2X2) / 8
  
  
  #B2C2
  B2C2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
  B2C2 <- -sum(xyCors2) * crossprod(Y2, crossprod(B2C2mat, Y3)) / 8
  
  
  #B2C3
  B2C3 <- -crossprod(Y2, crossprod(B2B2mat, Y * rX)) / 8
  
  
  
  #B2C4
  B2C4 <- -crossprod(Y2, crossprod(AB2mat, Y * r2X2)) / 8

  
  
  #B3^2
  B3B3 <- B2B2
  
  
  
  #B3B4
  B3B4 <- B2B4
  
  
  
  #B3C1
  B3C1 <- B2C2
  
  
  
  #B3C2
  B3C2 <- B2C1
  
  
  
  #B3C3
  B3C3 <- B2C4
  
  
  
  #B3C4
  B3C4 <- B2C3
  
  
  
  #B4^2
  B4B4 <- sum(B2B2mat^2) / 16
  
  
  
  #B4C1
  B4C1 <- -crossprod(r2X2, crossprod(B2C2mat, Y3)) / 8
  
  
  
  #B4C2
  B4C2 <- B4C1
  
  
  
  #B4C3
  B4C3 <- -sum(crossprod(Y, B2B2mat * B2C2mat)) / 8
  
  
  
  #B4C4
  B4C4 <- B4C3
  
  
  
  #C1^2
  C1C1 <- sum(xyCors2) * crossprod(crossprod(X, Y3)) / 4
  
  
  
  #C1C2
  C1C2 <- crossprod(Y3, rX)^2 / 4
  
  
  
  #C1C3
  C1C3 <- crossprod(crossprod(X, Y3), crossprod(X, Y * r2X2)) / 4
  
  
  
  #C1C4
  C1C4 <- crossprod(Y3, crossprod(t(AB2mat), Y * rX)) / 4
  
  
  
  #C2^2
  C2C2 <- C1C1
  
  
  
  #C2C3
  C2C3 <- C1C4
  
  
  
  #C2C4
  C2C4 <- C1C3
  
  
  
  #C3^2
  C3C3 <- crossprod(Y, crossprod(B2B2mat * XXt, Y)) / 4
  
  
  
  #C3C4
  C3C4 <- crossprod(Y, crossprod(AB2mat * t(AB2mat), Y)) / 4
  
  
  
  #C4^2
  C4C4 <- C3C3
  
  # Collect results
  resvec <- c(AA,  AB1,  AB2,  AB3,  AB4,  AC1,  AC2,  AC3,  AC4,
                  B1B1, B1B2, B1B3, B1B4, B1C1, B1C2, B1C3, B1C4,
                        B2B2, B2B3, B2B4, B2C1, B2C2, B2C3, B2C4, 
                              B3B3, B3B4, B3C1, B3C2, B3C3, B3C4,
                                    B4B4, B4C1, B4C2, B4C3, B4C4,
                                          C1C1, C1C2, C1C3, C1C4,
                                                C2C2, C2C3, C2C4, 
                                                      C3C3, C3C4,
                                                            C4C4)
  
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])
  
  ### Now calculating trace of the straight-up varmat
  star1 <- crossprod(Y2, X2RowSums)
  star2 <- Y4ColSums * sum(xyCors2)
  star3 <- crossprod(Y2, r2X2)
  star4 <- crossprod(xyCors2, X4ColSums)
  dagger1 <- crossprod(Y3, rX)
  dagger2 <- crossprod(Y, rX3)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2
  
  return(c(trmat / (n - 1), trprod / (n - 1)^2))
  
}

trace_fast <- function () {
  
  # General calcs
  xyCors <- as.vector(crossprod(Y, X)) / (n - 1)
  xyCors2 <- xyCors^2
  rX <- crossprod(tX, xyCors)
  r2X <- crossprod(tX, xyCors2)
  r2X2 <- crossprod(tX2, xyCors2)
  rX3 <- crossprod(tX3, xyCors)

  
  
  #AA
  AA <- crossprod(Y2, crossprod(XXt2, Y2))
  
  
  
  #AB1
  AB1 <- Y4ColSums * crossprod(Y2, rX^2) / 4
  
  
  
  #AB2
  AB2array <- do.call(cbind, lapply(1:n, function (i) tX[ , i] * tX2))
  AB2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
  AB2 <- crossprod(Y2, crossprod(AB2mat, Y2 * rX)) / 4
  
  #AB3
  AB3 <- AB2
  
  
  
  #AB4
  AB4 <- sum(Y2 * AB2mat^2) / 4
  
  
  
  #AC1
  AC1vec1 <- crossprod(Y2 * rX, X)
  AC1vec2 <- crossprod(X, Y3)
  AC1 <- -AC1vec1 %*% AC1vec2 / 2
  
  #AC2
  AC2 <- AC1
  
  
  
  #AC3
  AC3 <- -crossprod(Y, crossprod(AB2mat * XXt, Y2)) / 2
  
  
  #AC4
  AC4 <- AC3
  
  
  
  #B1^2
  B1B1 <- Y4ColSums^2 * sum(xyCors2)^2 / 16
  
  
  
  #B1B2
  B1B2 <- Y4ColSums * crossprod(Y2, r2X2) * sum(xyCors2) / 16
  
  
  
  #B1B3
  B1B3 <- B1B2
  
  
  
  #B1B4
  B1B4 <- Y4ColSums * crossprod(r2X2) / 16
  
  
  
  #B1C1
  B1C1 <- -Y4ColSums * sum(xyCors2) * crossprod(Y3, rX) / 8
  
  

  #B1C2
  B1C2 <- B1C1
  
  
  
  #B1C3
  B1C3 <- -Y4ColSums * crossprod(Y, rX * r2X2) / 8
  
  
  
  #B1C4
  B1C4 <- B1C3
  
  
  
  #B2^2
  B2B2array <- do.call(cbind, lapply(1:n, function (i) tX2[ , i] * tX2))
  B2B2mat <- matrix(crossprod(B2B2array, xyCors2), n, byrow = TRUE)
  B2B2 <- sum(xyCors2) * crossprod(Y2, crossprod(B2B2mat, Y2)) / 16
  
  #B2B3
  B2B3 <- crossprod(Y2, r2X2)^2 / 16
  
  
  
  #B2B4
  B2B4 <- crossprod(Y2, crossprod(B2B2mat, r2X2)) / 16
  
  
  
  #B2C1
  B2C1 <- -crossprod(Y3, rX) * crossprod(Y2, r2X2) / 8
  
  
  #B2C2
  B2C2mat <- matrix(crossprod(AB2array, xyCors), n, byrow = TRUE)
  B2C2 <- -sum(xyCors2) * crossprod(Y2, crossprod(B2C2mat, Y3)) / 8
  
  
  #B2C3
  B2C3 <- -crossprod(Y2, crossprod(B2B2mat, Y * rX)) / 8
  
  
  
  #B2C4
  B2C4 <- -crossprod(Y2, crossprod(AB2mat, Y * r2X2)) / 8

  
  
  #B3^2
  B3B3 <- B2B2
  
  
  
  #B3B4
  B3B4 <- B2B4
  
  
  
  #B3C1
  B3C1 <- B2C2
  
  
  
  #B3C2
  B3C2 <- B2C1
  
  
  
  #B3C3
  B3C3 <- B2C4
  
  
  
  #B3C4
  B3C4 <- B2C3
  
  
  
  #B4^2
  B4B4 <- sum(B2B2mat^2) / 16
  
  
  
  #B4C1
  B4C1 <- -crossprod(r2X2, crossprod(B2C2mat, Y3)) / 8
  
  
  
  #B4C2
  B4C2 <- B4C1
  
  
  
  #B4C3
  B4C3 <- -sum(crossprod(Y, B2B2mat * B2C2mat)) / 8
  
  
  
  #B4C4
  B4C4 <- B4C3
  
  
  
  #C1^2
  C1C1 <- sum(xyCors2) * crossprod(crossprod(X, Y3)) / 4
  
  
  
  #C1C2
  C1C2 <- crossprod(Y3, rX)^2 / 4
  
  
  
  #C1C3
  C1C3 <- crossprod(crossprod(X, Y3), crossprod(X, Y * r2X2)) / 4
  
  
  
  #C1C4
  C1C4 <- crossprod(Y3, crossprod(t(AB2mat), Y * rX)) / 4
  
  
  
  #C2^2
  C2C2 <- C1C1
  
  
  
  #C2C3
  C2C3 <- C1C4
  
  
  
  #C2C4
  C2C4 <- C1C3
  
  
  
  #C3^2
  C3C3 <- crossprod(Y, crossprod(B2B2mat * XXt, Y)) / 4
  
  
  
  #C3C4
  C3C4 <- crossprod(Y, crossprod(AB2mat * t(AB2mat), Y)) / 4
  
  
  
  #C4^2
  C4C4 <- C3C3
  
  # Collect results
  resvec <- c(AA,  AB1,  AB2,  AB3,  AB4,  AC1,  AC2,  AC3,  AC4,
                  B1B1, B1B2, B1B3, B1B4, B1C1, B1C2, B1C3, B1C4,
                        B2B2, B2B3, B2B4, B2C1, B2C2, B2C3, B2C4, 
                              B3B3, B3B4, B3C1, B3C2, B3C3, B3C4,
                                    B4B4, B4C1, B4C2, B4C3, B4C4,
                                          C1C1, C1C2, C1C3, C1C4,
                                                C2C2, C2C3, C2C4, 
                                                      C3C3, C3C4,
                                                            C4C4)
  
  resmat <- matrix(0, 9, 9)
  resmat[row(resmat) >= col(resmat)] <- resvec
  trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])
  
  ### Now calculating trace of the straight-up varmat
  star1 <- crossprod(Y2, X2RowSums)
  star2 <- Y4ColSums * sum(xyCors2)
  star3 <- crossprod(Y2, r2X2)
  star4 <- crossprod(xyCors2, X4ColSums)
  dagger1 <- crossprod(Y3, rX)
  dagger2 <- crossprod(Y, rX3)
  trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2
  
  return(c(trmat / (n - 1), trprod / (n - 1)^2))
  
}

