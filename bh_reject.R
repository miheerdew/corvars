library(Ckmeans.1d.dp)

symdiff <- function (s1, s2) {
  return(union(setdiff(s1, s2), setdiff(s2, s1)))
}

jaccard <- function (s1, s2) {
  return(length(symdiff(s1, s2)) / length(union(s1, s2)))
}

cluster_thres <- function (zs) {
  
  clust <- Ckmeans.1d.dp(zs, 2)$cluster
  thres <- min(zs[clust == 2])
  return(which(zs >= thres))
}

bh_reject <- function (pvals, alpha, conserv = TRUE) {
  
  m <- length(pvals)
  
  if (!conserv) {
    pvals_adj <- m * pvals / rank(pvals)
  } else {
    mults <- sum(1 / c(1:m))
    pvals_adj <- mults * m * pvals / rank(pvals)
  }
  
  if (sum(pvals_adj <= alpha) > 0) {
    thres <- max(pvals[pvals_adj <= alpha])
    return(which(pvals <= thres))
  } else {
    return(integer(0))
  }
}

bhy <-
  function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bhy cutoff
    nums = 1:length(pvals)
    cms = sum(1/nums)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals)*cms)*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
      return(c())
    }else{
      cutoff = max(which(under))
      return(ord[1:cutoff])
    }
  }

bh <-
  function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bh cutoff
    nums = 1:length(pvals)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals))*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
      return(c())
    }else{
      cutoff = max(which(under))
      return(ord[1:cutoff])
    }
  }
