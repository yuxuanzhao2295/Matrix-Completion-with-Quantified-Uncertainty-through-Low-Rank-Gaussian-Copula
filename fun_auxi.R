generate_LRGC = function(n=500,p.seq=c(100,100,100),rank,sigma, seed = NULL,
                         ordinal.new = FALSE, ord.num = 5,
                         cont.type = 'LR',
                         fun = NULL, Wreturn = FALSE){
  if (!is.null(seed)) set.seed(seed)
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
      if (ordinal.new) X_true[,m] = continuous2ordinal_new(Z[,m], k = 2)
      else X_true[,m] = continuous2ordinal(Z[,m], k = 2)
    }
  }
  if (!is.null(index.o)){
    for (m in index.o){
      if (ordinal.new) X_true[,m] = continuous2ordinal_new(Z[,m], k = ord.num)
      else X_true[,m] = continuous2ordinal(Z[,m], k = ord.num)
    }
  }
  if (Wreturn) list(X = X_true, W = W) else X_true
}