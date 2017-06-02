library(emulator)

pvals_chisq <- function (S) {
  
  testY <- (S[1]) <= dx
  
  # General setup
  
  if (testY) {
    
    B <- S
    X <- X0$value[ , B]
    Y <- copy(Y0)
    
    X2 <- X02$value[ , B]
    X3 <- X03$value[ , B]
    X4 <- X04$value[ , B]
    Y2 <- copy(Y02)
    Y3 <- copy(Y03)
    Y4 <- copy(Y04)
    
    tX <- tX0$value[B, ]
    tY <- copy(tY0)
    tX2 <- tX02$value[B, ]
    tX3 <- tX03$value[B, ]
    
    Y4ColSums <- copy(Y04ColSums)
    X4ColSums <- X04ColSums$value[B]
    
    allr <- copy(xyCors)
    allr22 <- copy(xyCors22)
    allr31 <- copy(xyCors31_y)
    
  } else {
    
    X <- Y0$value[ , B]
    Y <- copy(X0)
    
    X2 <- Y02$value[ , B]
    X3 <- Y03$value[ , B]
    X4 <- Y04$value[ , B]
    Y2 <- copy(X02)
    Y3 <- copy(X03)
    Y4 <- copy(X04)
    
    tX <- tY0$value[B, ]
    tY <- copy(tX0)
    tX2 <- tY02$value[B, ]
    tX3 <- tY03$value[B, ]
    
    Y4ColSums <- copy(X04ColSums)
    X4ColSums <- Y04ColSums$value[B]
    
    allr <- copy(txyCors)
    allr22 <- copy(txyCors22)
    allr31 <- copy(xyCors31_x)
    
  }
  
  
  if (length(B) > n) { # Do n^2|B| calc
    
    XXt <- tcrossprod(X)
    XXt2 <- XXt^2
    X1RowSums <- rowSums(X)
    X2RowSums <- rowSums(X2)
    
    trs <- sapply(1:ncol(Y$value), function (i) {
        
      m <- n
      y <- Y$value[ , i]
      e <- matrix(1, nrow=m, ncol=1)
      y2 <- Y2$value[ , i] # m x 1  ## GET THIS FROM Y2
      y3 <- Y3$value[ , i] ## GET THIS FROM Y3
      
      M <- m - 1
      r <- allr$value[i, B, drop = FALSE] / M # 1 x dx ## GET THIS FROM allr
      r2 <- allr22$value[i, B, drop = FALSE] / M # 1 x dx
      r3 <- allr31$value[i, B, drop = FALSE] / M # 1 x dx
      rX <- tcrossprod(X, r) # m x 1
      r3X <- tcrossprod(X, r3) # m x 1
      rX3 <- tcrossprod(X3, r) # m x 1
      
      
      rSq <- r^2
      rSqX2 <- tcrossprod(X2,rSq) # m x 1
      r2rSqX2 <- tcrossprod(X2, r2 * rSq) # m x 1
      r3rX2 <- tcrossprod(X2, r3*r) # m x 1
      
      tXr <- tX * as.vector(r) # dx x m
      
      y4Sum <- Y4ColSums$value[i] # scalar ## GET THIS FROM Y4ColSums
      rSqSum <- sum(rSq) #scalar
      
      Z1 <- X2 %*% tXr # m x m
      Z2 <- X2 %*% tXr^2 # m x m
      
      #All of these are scalars
      #TODO: Normalize them by (m-1)^2 after they are working correctly.
      
      AA <- quad.form(XXt2, y2)
      
      AB1 <- y4Sum * crossprod(y2, rX^2)
      
      AB2 <- quad.3form(XXt, y, y2 * tcrossprod(X, r2*r))
      
      AB3 <- AB2
      
      AB4 <- quad.3form(Z1^2, e, y2)
      
      AC1 <- crossprod(rX * r3X, y2) * M
      
      AC2 <- quad.3form(Z1 * XXt, y, y2)
      
      AC3 <- AC1
      
      AC4 <- AC2
      
      B1B1 <- (y4Sum*rSqSum)^2
      
      B1B2 <- rSqSum * y4Sum * tcrossprod(r2,rSq) * M
      
      B1B3 <- B1B2
      
      B1B4 <- y4Sum * sum(rSqX2^2)
      
      B1C1 <- y4Sum * rSqSum * tcrossprod(r,r3) * M
      
      B1C2 <- y4Sum * crossprod(y, rX * rSqX2)
      
      B1C3 <- B1C1
      
      B1C4 <- B1C2
      
      B2B2 <- tcrossprod(rSq, r2^2) * rSqSum * M^2
      
      B2B3 <- tcrossprod(rSq, r2)^2 * M^2
      
      B2B4 <- sum(rSqX2 * r2rSqX2) * M
      
      B2C1 <- tcrossprod(rSq, r2) * tcrossprod(r,r3) * M^2
      
      B2C2 <- crossprod(y, r2rSqX2 * rX) * M
      
      B2C3 <- sum(r * r2 * r3) * rSqSum * M^2
      
      B2C4 <- crossprod(y, tcrossprod(X, r * r2) * rSqX2) * M
      
      B3B3 <- B2B2
      
      B3B4 <- B2B4
      
      B3C1 <- B2C3
      
      B3C2 <- B2C4
      
      B3C3 <- B2C1
      
      B3C4 <- B2C2
      
      
      
      B4B4 <- sum(Z2^2)
      
      B4C1 <- sum(rSqX2 * r3rX2) * M
      
      B4C2 <- quad.3form(Z2 * Z1, e, y)
      
      B4C3 <- B4C1
      
      B4C4 <- B4C2
      
      C1C1 <- rSqSum * sum(r3^2) * M^2
      
      C1C2 <- crossprod(y, rSqX2 * r3X) * M
      
      C1C3 <- sum(r * r3)^2 * M^2
      
      C1C4 <- crossprod(y, rX * r3rX2) * M
      
      C2C2 <- quad.form(Z2 * XXt, y)
      
      C2C3 <- C1C4
      
      C2C4 <- quad.form(Z1*t(Z1), y)
      
      C3C3 <- C1C1
      
      C3C4 <- C1C2
      
      C4C4 <- C2C2
      
      resvec <- c(AA, c(AB1, AB2,  AB3,  AB4)/4, -1/2*c( AC1,  AC2,  AC3,  AC4),
                  c(B1B1, B1B2, B1B3, B1B4)/16, -1/8*c(B1C1, B1C2, B1C3, B1C4),
                  c(B2B2, B2B3, B2B4)/16, -1/8*c(B2C1, B2C2, B2C3, B2C4),
                  c(B3B3, B3B4)/16, -1/8*c(B3C1, B3C2, B3C3, B3C4),
                  B4B4/16, -1/8*c(B4C1, B4C2, B4C3, B4C4),
                  c(C1C1, C1C2, C1C3, C1C4,
                    C2C2, C2C3, C2C4,
                    C3C3, C3C4,
                    C4C4)/4)
      resmat <- matrix(0, 9, 9)
      resmat[row(resmat) >= col(resmat)] <- resvec
      trprod <- sum(diag(resmat)) + 2 * sum(resmat[row(resmat) > col(resmat)])
      
      #Trace of straight up varmat
      star1 <- crossprod(y2, X2RowSums) #O(m)
      star2 <- y4Sum * rSqSum 
      star3 <- crossprod(y2, rSqX2) #O(m) 
      star4 <- tcrossprod(rSq, X4ColSums) #O(n)
      dagger1 <- crossprod(y3, rX) #O(m)
      dagger2 <- crossprod(y, rX3) #O(m)
      trmat <- star1 + (star2 + 2 * star3 + star4) / 4 - dagger1 - dagger2
      
      return(c(trmat/M, trprod/M^2))
      
    })
    
    
  } else {            # Do n|B|^2 calc
    
    # General calcs
    xyCors <- as.vector(allr$value[i, B, drop = FALSE]) / M
    rX2 <- crossprod(tX2, xyCors)
    
    # star 1
    S1 <- crossprod(as.vector(Y$value[ , i]) * X)
    
    # star 2
    S2 <- Y4ColSums$value[i] * tcrossprod(xyCors)
    
    # star 3
    S3 <- tcrossprod(xyCors, xyCors * colSums(X2 * as.vector(Y2$value[ , i])))
    
    # star 4
    S4 <- tcrossprod(xyCors * t(X2))
    
    # dagger 1
    D1 <- tcrossprod(xyCors, colSums(X * as.vector(Y3$value[ , i])))
    
    # dagger 2
    D2 <- t(crossprod(as.vector(Y$value[ , i]) * X, X2)) * xyCors
    
    bigmat <- (S1 + (S2 + S3 + t(S3) + S4) / 4 - (D1 + t(D1) + D2 + t(D2)) / 2) / (n - 1)
    
    trSig <- sum(diag(bigmat))
    trSig2 <- sum(diag(crossprod(bigmat)))
    
    return(c(trSig, trSig2))
    
  }
  
}