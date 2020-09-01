ppca_edited <- function(Matrix, nPcs=2, seed=NA, threshold=1e-5, maxIterations=1000, ...) {
  ## Set the seed to the user defined value. This affects the generation
  ## of random values for the initial setup of the loading matrix
  if (!is.na(seed)) 
    set.seed(seed)
  
  N <- nrow(Matrix)
  D <- ncol(Matrix)
  
  Obs <- !is.na(Matrix)
  hidden <- which(is.na(Matrix))
  missing <- length(hidden)
  
  if(missing) { Matrix[hidden] <- 0 } 
  
  ## ------- Initialization
  r <- sample(N)
  C <- t(Matrix[r[1:nPcs], ,drop = FALSE])
  ## Random matrix with the same dimnames as Matrix
  C <- matrix(rnorm(C), nrow(C), ncol(C), dimnames = labels(C) )
  CtC <- t(C) %*% C
  ## inv(C'C) C' X is the solution to the EM problem
  X  <- Matrix %*% C %*% solve(CtC)
  recon    <- X %*% t(C)
  recon[hidden] <- 0
  ss <- sum(sum((recon - Matrix)^2)) / (N * D - missing)
  
  count <- 1
  old <- Inf
  
  ## ------ EM iterations
  while (count > 0) {
    ## E-step, (co)variances
    Sx <- solve(diag(nPcs) + CtC/ss) 
    ss_old <- ss
    if(missing) {
      proj <- X %*% t(C)
      Matrix[hidden] <- proj[hidden]
    }
    
    ## E step: expected values
    X <- Matrix %*% C %*% Sx / ss
    
    ## M-step
    SumXtX <- t(X) %*% X
    
    ## Replace the right matrix division from matlab
    C <- (t(Matrix) %*% X) %*% solve( (SumXtX + N * Sx) )
    
    CtC <- t(C) %*% C
    ss <- ( sum(sum( (C %*% t(X) - t(Matrix))^2 )) + N * sum(sum(CtC %*% Sx)) +
              missing * ss_old ) / (N * D)
    
    objective <- N * (D * log(ss) + sum(diag(Sx)) - log(det(Sx)) ) +
      sum(diag(SumXtX)) - missing * log(ss_old)
    
    rel_ch <- abs( 1 - objective / old )
    old <- objective
    
    count <- count + 1
    if( rel_ch < threshold & count > 5 ) {
      count <- 0
    }
    else if (count > maxIterations) {
      count <- 0
      warning("stopped after max iterations, but rel_ch was > threshold")
    }
  } ## End EM iteration
  C <- orth(C)
  evs <- eigen( cov(Matrix %*% C) )
  vals <- evs[[1]]
  vecs <- evs[[2]]
  
  C <- C %*% vecs
  X <- Matrix %*% C
  
  ## Paramters in original Matlab implementation were:
  ## C (D by d)    - C has the approximate loadings (eigenvectors of
  ## the covariance matrix)
  ##          as columns.
  ## X        - The approximate scores 
  ## Matrix (N by D)    - Expected complete observations.
  ## M (D by 1)    - Column wise data mean
  ## ss (scalar)    - isotropic variance outside subspace
  
  R2cum <- rep(NA, nPcs)
  TSS <- sum(Matrix^2, na.rm=TRUE)
  for (i in 1:ncol(C)) {
    difference <- Matrix - (X[,1:i, drop=FALSE] %*% t(C[,1:i, drop=FALSE]))
    R2cum[i] <- 1 - (sum(difference^2, na.rm=TRUE) / TSS)
  }
  
  res <- new("pcaRes")
  res@scores <- X
  res@loadings <- C
  res@R2cum <- R2cum
  res@method <- "ppca"
  
  cObs <- Matrix
  cObs[missing] <- fitted(res, post=TRUE)[missing]
  res@completeObs <- cObs
  return(list(res = res, sigma = ss))
}

orth <- function(mat, skipInac = FALSE) {
  
  if(nrow(mat) > ncol(mat)) {
    leftSVs <- ncol(mat)
  } else {
    leftSVs <- nrow(mat)
  }
  
  result <- svd(mat, nu =  leftSVs, nv = ncol(mat))
  U <- result[[2]]
  S <- result[[1]]
  V <- result[[3]]
  
  m <- nrow(mat)
  n <- ncol(mat)
  
  if(m > 1) { 
    s <- diag(S, nrow = length(S))
  } else     if(m == 1) { 
    s <- S[1] 
  } else { 
    s <- 0 
  }
  
  tol <- max(m,n) * max(s) * .Machine$double.eps
  r <- sum(s > tol)
  if ( r < ncol(U) ) {
    if (skipInac) {
      warning("Precision for components ", r + 1 , " - ", ncol(U), 
              " is below .Machine$double.eps. \n",
              "Results for those components are likely to be inaccurate!!\n",
              "These component(s) are not included in the returned solution!!\n")
    } else {
      warning("Precision for components ", r + 1 , " - ", ncol(U), 
              " is below .Machine$double.eps. \n",
              "Results for those components are likely to be inaccurate!!\n")
    }
  }
  
  if (skipInac) {
    ONB <- U[, 1:r, drop=FALSE]
    ## Assing correct row and colnames
    rownames(ONB) <- labels(mat[, 1:r, drop=FALSE])[[1]];
    colnames(ONB) <- labels(mat[, 1:r, drop=FALSE])[[2]];
  } else {
    ONB<-U
    ## Assing correct row and colnames
    rownames(ONB) <- labels(mat)[[1]];
    colnames(ONB) <- labels(mat)[[2]];
  }
  
  return(ONB)
}

ct_impute_ppca = function(X, fit, alpha = 0.95){
  # output: confidence interval for predicted continuous entries
  # fit is a returned list from function impute_mixedgc_ppca
  # used objects: W, sigma, Zimp (imputed values of Z)
  # The original observation X is used to detect the missing locations and estimate the marginals 
  n = dim(X)[1]
  p = dim(X)[2]
  obj = svd(fit$W)
  U = obj$u
  d = obj$d
  rank = length(d)
  
  sig = fit$sigma
  Zimp = fit$Zimp
  
  upper = array(NA, dim = c(n,p))
  lower = array(NA, dim = c(n,p))
  var = array(NA, dim = c(n,p))
  margin = qnorm(1-(1-alpha)/2)
  # compute conditional variance 
  for (i in 1:n){
    index_m = which(is.na(X[i,]))
    index_o = which(!is.na(X[i,]))
    Uiobs = matrix(U[index_o,], ncol = rank)
    Uimis = matrix(U[index_m,], ncol = rank)
    
    dUmis = solve(sig * diag(d^{-2}) + t(Uiobs) %*% Uiobs, t(Uimis))
    # compute variance i.e. diagonal elements of conditional covariance matrix
    for (l in 1:length(index_m)){
      j = index_m[l]
      du = dUmis[,l]
      var_ij = sig + sig * sum(du * U[j,]) 
      var[i,j] = var_ij
      upper[i,j] = Zimp[i,j] + margin * sqrt(var_ij)
      lower[i,j] = Zimp[i,j] - margin * sqrt(var_ij)
    }
  }
  
  
  return(list(upper = upper, lower = lower, var = var))
}