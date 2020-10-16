generate_LRGC = function(n=500,p.seq=c(100,100,100),rank,sigma, seed = 1,
                         ord.num = 5,
                         cont.type = 'LR',
                         fun = NULL, Wreturn = FALSE){
  set.seed(seed)
  require(mvtnorm)
  require(mixedgcImp)
  #
  index.c = NULL
  index.b = NULL
  index.o = NULL
  if (p.seq[1]>0) index.c = seq(1,p.seq[1],by = 1)
  if (p.seq[2]>0) index.o = seq(1,p.seq[2],by = 1) + p.seq[1]
  if (p.seq[3]>0) index.b = seq(1,p.seq[3],by = 1) + p.seq[1] + p.seq[2]
  p = sum(p.seq)
  
  W = matrix(rnorm(p*rank), nrow = p)
  W = t(apply(W, 1, function(x){x/sqrt(sum(x^2))  * sqrt((1-sigma))}))
  Z = matrix(rnorm(n*rank), ncol = rank) %*% t(W) + matrix(rnorm(n*p, sd = sqrt(sigma)), ncol = p)
  X_true = Z
  #
  if (!is.null(index.c)){
    if (is.null(fun)){
      if (cont.type == 'LR')  fun = function(x){x} else fun = function(x){x^3}
    }
    for (m in index.c) X_true[,m] = fun(Z[,m])
  }
  if (!is.null(index.b)){
    for (m in index.b){
      X_true[,m] = continuous2ordinal_new(Z[,m], k = 2)
    }
  }
  if (!is.null(index.o)){
    for (m in index.o){
      X_true[,m] = continuous2ordinal_new(Z[,m], k = ord.num)
    }
  }
  if (Wreturn) list(X = X_true, W = W) else X_true
}


data_split = function(X, ratio.val = 0.1, ratio.test = 0.1, seed = 101){
  loc = which(!is.na(X))
  i= 0
  empty_row = TRUE
  while(empty_row){
    set.seed(seed)
    loc_val = sample(loc, round(length(loc) * ratio.val))
    loc_test = sample(setdiff(loc,loc_val), round(length(loc) * ratio.test))
    loc_train = setdiff(loc, union(loc_test,loc_val))
    
    Xtrain = X
    Xtrain[union(loc_test,loc_val)] = NA
    
    empty_row = any(apply(Xtrain, 1, function(x){sum(!is.na(x))}) == 0)
    i=i+1
    seed=seed+1
    if (i > 100) stop('cannot produce masking without empty row')
  }
  return(list(Xtrain = Xtrain, seed=seed-1, 
              loc_test = loc_test, loc_train = loc_train, loc_val = loc_val))
}

continuous2ordinal_new = function (x, k = 2, cutoff = NULL){
  q = quantile(x, c(0.05,0.95))
  if (k == 2) {
    if (is.null(cutoff)) {
      cutoff = sample(x[x>q[1] & x<q[2]],1)
    }
    x = (x >= cutoff)
  }
  else {
    if (is.null(cutoff)) {
      #cutoff = seq(min(x) - 0.1 * sd(x), max(x) + 0.1 *  sd(x), length.out = k + 1)
      cutoff = c(min(x) - 0.1 * sd(x), sort(sample(x[x>q[1] & x<q[2]], k-1)),max(x) + 0.1 * sd(x))
    }
    x = cut(x, cutoff, labels = FALSE, include.lowest = TRUE)
  }
  x
}

mae_by_reliability = function(r,xtrue,ximp,l=100,start=1){
  # r for reliability
  q = seq(0,1-start/l,by = start/l)
  r_q = quantile(r, q)
  
  err = numeric(l)
  for (i in 1:l){
    loc_q = which(r >= r_q[i])
    
    val = xtrue[loc_q]
    imp = ximp[loc_q]
    
    err[i] = mean(abs(val - imp)) 
  }
  
  err
}

nrmse_by_reliability = function(r,xtrue,ximp,l=100,start=1){
  # r for reliability
  q = seq(0,1-start/l,by = start/l)
  r_q = quantile(r, q)
  err = numeric(l)
  for (i in 1:l){
    loc_q = which(r >= r_q[i])
    
    val = xtrue[loc_q]
    imp = ximp[loc_q]
    
    err[i] = sqrt(mean((val - imp)^2)/mean(val^2))
  }
  
  err
}

rmse_by_reliability = function(r,xtrue,ximp,l=100,start=1){
  q = seq(0,1-start/l,by = start/l)
  r_q = quantile(r, q)
  err = numeric(l)
  for (i in 1:l){
    loc_q = which(r >= r_q[i])
    
    val = xtrue[loc_q]
    imp = ximp[loc_q]
    
    err[i] = sqrt(mean((val - imp)^2))
  }
  
  err
}

trunc.rating = function(x, xmin = 1, xmax = 5){
  x = round(x)
  x[which(x<xmin)] = xmin
  x[which(x>xmax)] = xmax
  x
}

varest_softImpute = function(X, loc = NULL, ratio, nfolds = 10, seed = 1, parallel = FALSE, maxrank = dim(X)[2]-1){
  # mask additional and fit 10 more times
  if (is.null(loc)) loc = which(is.na(X))
  
  m = CV_split(X, nfolds = nfolds, seed = seed)
  if (parallel){
    require(doParallel)
    imp_cv = foreach (s = 1:nfolds, .combine = rbind)%dopar%{
      Xnew = X
      Xnew[m[[s]]] = NA # additional masking
      xc = biScale(Xnew, row.scale=FALSE, col.scale=FALSE)
      lam0 = lambda0(xc) * ratio
      fiti = softImpute(xc, lambda=lam0, rank=maxrank, type = 'svd')
      complete(xc, fiti)[loc]
    }
  }else{
    imp_cv = array(0,dim = c(length(loc),nfolds)) 
    for (s in 1:nfolds){
      Xnew = X
      Xnew[m[[s]]] = NA # additional masking
      xc = biScale(Xnew, row.scale=FALSE, col.scale=FALSE)
      lam0 = lambda0(xc) * ratio
      fiti = softImpute(xc, lambda=lam0, rank=maxrank, type = 'svd')
      imp_cv[,s] = complete(xc, fiti)[loc]
    }
  }
  
  imp_var = apply(imp_cv, 1, var)
  imp_var
}

CV_split = function(xdata, nfolds = 10, seed = 1){
  # generate location
  ind = which(!is.na(xdata))
  size = round(length(ind)/nfolds)
  set.seed(seed)
  ind = sample(ind, length(ind))
  loc.list = list()
  for (i in 1:nfolds){
    start = (i-1)*size + 1
    if (i<nfolds) end = i*size else end = length(ind)
    loc.list[[i]] = ind[start:end]
  }
  loc.list
}