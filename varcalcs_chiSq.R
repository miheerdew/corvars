varcalc1 <- function (Yvec, Xmat) {
  
  Xmat <- scale(Xmat)
  Yvec <- as.vector(scale(Yvec))
  
  n <- nrow(Xmat)
  
  if (n != length(Yvec))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- as.vector(cor(Yvec, Xmat))
  
  y4 <- sum(Yvec^4)
  x2RowSums <- rowSums(Xmat^2)
  x4ColSums <- colSums(Xmat^4)
  
  corX <- tcrossprod(xyCors, Xmat)
  cor2X2 <- tcrossprod(xyCors^2, Xmat^2)
  corX3 <- tcrossprod(xyCors, Xmat^3)
  
  # Calc for star 1
  star1 <- crossprod(Yvec^2, x2RowSums)
  
  # Calc for star 2
  star2 <- y4 * sum(xyCors^2)
  
  # Calc for star 3
  star3 <- 2 * crossprod(Yvec^2, t(cor2X2))
  
  # Calc for star 4
  star4 <- crossprod(xyCors^2, x4ColSums)
  
  # Calc for dagger 1
  dagger1 <- tcrossprod(Yvec^3, corX)
  
  # Calc for dagger 2
  dagger2 <- tcrossprod(Yvec, corX3)
  
  
  
  return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}

varcalc1_multi <- function (Ymat, Xmat) {
  
  Xmat <- scale(Xmat)
  Ymat <- scale(Ymat)
  
  n <- nrow(Xmat)
  
  if (n != nrow(Ymat))
    stop('X and Y variables must be of same dimension\n')
  
  # General calcs
  xyCors <- cor(Ymat, Xmat)
  
  y4 <- colSums(Ymat^4)
  x2RowSums <- rowSums(Xmat^2)
  x4ColSums <- colSums(Xmat^4)
  
  corX <- tcrossprod(xyCors, Xmat)
  cor2X2 <- tcrossprod(xyCors^2, Xmat^2)
  corX3 <- tcrossprod(xyCors, Xmat^3)
  
  # Calc for star 1
  star1 <- crossprod(Ymat^2, x2RowSums)
  
  # Calc for star 2
  star2 <- y4 * rowSums(xyCors^2)
  
  # Calc for star 3
  star3 <- 2 * rowSums(t(Ymat^2) * cor2X2)
  
  # Calc for star 4
  star4 <- crossprod(t(xyCors^2), x4ColSums)
  
  # Calc for dagger 1
  dagger1 <- rowSums(t(Ymat^3) * corX)
  
  # Calc for dagger 2
  dagger2 <- rowSums(t(Ymat) * corX3)
  
    return((star1 + 0.25 * (star2 + star3 + star4) - dagger1 - dagger2) / (n - 1))
  
  
}