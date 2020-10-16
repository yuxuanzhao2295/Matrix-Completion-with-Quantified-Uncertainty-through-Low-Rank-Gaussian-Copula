using LowRankModels;
using DataFrames;
using Random;
using Statistics;
using Printf;
using RCall;
using CSV;
R"library(mvtnorm)";
R"library(mixedgcImp)";
R"library(softImpute)";

# CALLED R FUNCTIONS 
R_generate_LRGC = R"generate_LRGC = function(n=500,p.seq=c(100,100,100),rank,sigma, seed = 1,
                         ord.num = 5,cont.type = 'LR',fun = NULL){
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
 
  X_true
}";

R"continuous2ordinal_new = function (x, k = 2, cutoff = NULL){
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
}"

R_mask = R"mask = function(Xnew, ratio, seed = 1){
  set.seed(seed)
  n = dim(Xnew)[1]
  p = dim(Xnew)[2]
  l = sample(1:(n*p), n*p*ratio)
  Xnew[l] = NA
  loc_mis = arrayInd(l, .dim = c(n,p))
  loc_obs = arrayInd(setdiff(1:(n*p),l), .dim=c(n,p))
  return(list(loc_mis = loc_mis, loc_obs = loc_obs, Xobs = Xnew))
}";

R_CV_split = R"CV_split = function(xdata,  nfold = 10, seed = 1){
  # generate location
  n = dim(xdata)[1]
  p = dim(xdata)[2]
  ind = which(!is.na(xdata))
  size = round(length(ind)/nfold)
  set.seed(seed)
  ind = sample(ind, length(ind))
  obs.1d = list()
  mis.1d = list()
  for (i in 1:nfold){
    start = (i-1)*size + 1
    if (i<nfold) end = i*size else end = length(ind)
    obs.1d[[i]] = setdiff(ind, ind[start:end])
    mis.1d[[i]] = setdiff(1:(n*p),obs.1d[[i]])
  }
  obs.2d = list()
  mis.2d = list()
  for (i in 1:nfold){
    obs.2d[[i]] = arrayInd(obs.1d[[i]], .dim = c(n,p))
    mis.2d[[i]] = arrayInd(mis.1d[[i]], .dim = c(n,p))
  }
  return(list(obs.1d = obs.1d, mis.1d = mis.1d, obs.2d = obs.2d, mis.2d = mis.2d))
}";

R_mae_by_reliability = R"mae_by_reliability = function(r,xtrue,ximp,l=100,start=1){
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
};"

R_nrmse_by_reliability = R"nrmse_by_reliability = function(r,xtrue,ximp,l=100,start=1){
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
}"

R_rmse_by_reliability = R"rmse_by_reliability = function(r,xtrue,ximp,l=100,start=1){
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
}"



R_data_split = R"data_split = function(X, ratio.val = 0.1, ratio.test = 0.1, seed = 101){
  n = dim(X)[1]
  p = dim(X)[2]
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
  loc_test_1d = loc_test
  loc_train_1d = loc_train
  loc_val_1d = loc_val
  loc_test = arrayInd(loc_test, .dim = c(n,p))
  loc_train = arrayInd(loc_train, .dim = c(n,p))
  loc_val = arrayInd(loc_val, .dim = c(n,p))
  return(list(loc_test = loc_test, loc_train = loc_train, loc_val = loc_val, Xtrain = Xtrain,
              loc_test_1d = loc_test_1d, loc_train_1d = loc_train_1d, loc_val_1d = loc_val_1d))
}";

# USED GLRM FUNCTIONS
function imp_MMMF_BvS(lambda, loc_obs, loc_mis, data)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    k = 199;
    @timed begin
        glrm = GLRM(data, BvSLoss(5), rx, ry, k, obs = loc_obs, scale = true);
        init_svd!(glrm);
        fit!(glrm, ProxGradParams(max_iter = 2000));
        imp = impute_missing(glrm);
        Val_test = [data[i[1], i[2]] for i in loc_mis];
        Val_imp = [imp[i[1], i[2]] for i in loc_mis];
        mean(abs.(Val_imp - Val_test))
    end
end


function imp_MMMF_l1(lambda, loc_obs, loc_mis, data)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    k =199;
    @timed begin
        glrm = GLRM(data, L1Loss(), rx, ry, k, obs = loc_obs, scale = true);
        init_svd!(glrm);
        fit!(glrm, ProxGradParams(max_iter = 2000));
        imp = impute_missing(glrm);
        imp = [max(min(round(i),5),1) for i in imp];
        Val_test = [data[i[1], i[2]] for i in loc_mis];
        Val_imp = [imp[i[1], i[2]] for i in loc_mis];
        mean(abs.(Val_imp - Val_test))
    end
end

function imp_MMMF_l2(lambda, loc_obs, loc_mis, data)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    k =199;
    @timed begin
        glrm = GLRM(data, QuadLoss(), rx, ry, k, obs = loc_obs, scale = true);
        init_svd!(glrm);
        fit!(glrm, ProxGradParams(max_iter = 2000));
        imp = impute_missing(glrm);
        Val_test = [data[i[1], i[2]] for i in loc_mis];
        Val_imp = [imp[i[1], i[2]] for i in loc_mis];
        mean((Val_imp - Val_test).^2)^.5/mean(Val_test.^2)^.5
    end
end

function imp_MMMF_binary(lambda, loc_obs, loc_mis, data, loss)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    k =199;
    @timed begin
        glrm = GLRM(data, loss, rx, ry, k, obs = loc_obs, scale = true);
        init_svd!(glrm);
        fit!(glrm, ProxGradParams(max_iter = 2000));
        imp = impute_missing(glrm);
        Val_test = [data[i[1], i[2]] for i in loc_mis];
        Val_imp = [imp[i[1], i[2]] for i in loc_mis];
        Val_test = convert(Array{Int8}, Val_test);
        Val_imp = convert(Array{Int8}, Val_imp);
        mean(abs.(Val_imp - Val_test))
    end
end

function imp_MMMF(loss, lambda, loc_obs, loc_mis, data, toint=false, k=199, maxiter = 2000)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    
    glrm = GLRM(data, loss, rx, ry, k, obs = loc_obs, scale = true);
    init_svd!(glrm);
    fit!(glrm, ProxGradParams(max_iter = 2000));
    imp = impute_missing(glrm);
    
    Val_test = [data[i[1], i[2]] for i in loc_mis];
    Val_imp = [imp[i[1], i[2]] for i in loc_mis];
    
    if toint
        Val_test = convert(Array{Int8}, Val_test);
        Val_imp = convert(Array{Int8}, Val_imp);
    end
    
    
    Val = zeros(length(Val_imp), 2);
    Val[:,1] = Val_imp;
    Val[:,2] = Val_test;
    Val
end

function imp_MMMF_XY(loss, lambda, loc_obs, loc_mis, data, k=199, maxiter = 2000)
    l_obs, l_mis = size(loc_obs)[1], size(loc_mis)[1];
    loc_obs = [(loc_obs[j,1], loc_obs[j,2]) for j in 1:l_obs];
    loc_mis = [(loc_mis[j,1], loc_mis[j,2]) for j in 1:l_mis];
    
    ry = QuadReg(lambda);
    rx = QuadReg(lambda);
    
    glrm = GLRM(data, loss, rx, ry, k, obs = loc_obs, scale = true);
    init_svd!(glrm);
    fit!(glrm, ProxGradParams(max_iter = maxiter));
    
    val = transpose(glrm.X) * glrm.Y;
    val = [val[i[1],i[2]] for i in loc_mis]
    val
end






