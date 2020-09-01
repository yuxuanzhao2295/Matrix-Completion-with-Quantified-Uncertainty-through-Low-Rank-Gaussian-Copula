ct_LRMC = function(u, d, v, X_obs){
  r = sum(round(d,4)>0)
  n = dim(X_obs)[1]
  p = dim(X_obs)[2]
  loc = which(is.na(X_obs))
  Ximp = u %*% diag(d) %*% t(v)
  # estimate noise
  sig = mean((Ximp[-loc] - X_obs[-loc])^2)
  # debiase
  Ximp[-loc] = Ximp[-loc] - (Ximp[-loc] - X_obs[-loc]) / ratio
  # projection 
  svd = svd(Ximp)
  X = svd$u[,1:r] %*% diag(sqrt(svd$d[1:r]))
  Y = svd$v[,1:r] %*% diag(sqrt(svd$d[1:r]))
  XX = solve(t(X) %*% X, t(X))
  YY = solve(t(Y) %*% Y, t(Y))
  
  upper = array(NA, dim = c(n,p))
  lower = array(NA, dim = c(n,p))
  
  len = numeric(length(loc))
  cov = numeric(length(loc))
  for (s in 1:length(loc)){
    ss = loc[s]
    coor = arrayInd(ss,.dim = c(n,p))
    i1 = coor[1]
    j1 = coor[2]
    vij = (sum(X[i1,] * XX[,i1]) + sum(Y[j1,] * YY[,j1])) * sig / ratio
    margin = sqrt(vij) * qnorm(1-0.025)
    upper[ss] = Ximp[ss] + margin
    lower[ss] = Ximp[ss] - margin
  }
  
  list(upper = upper, lower = lower)
}