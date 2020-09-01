# simulation in producing confidence interval
library(mixedgcImp)
library(softImpute)  # required for LRMC method, replicates the method in Chen, Yuxin, et al. "Inference and uncertainty quantification for noisy matrix completion." Proceedings of the National Academy of Sciences 116.46 (2019): 22931-22937.
library(missMDA) # required for MIPCA method
library(pcaMethods) # required for PPCA method
# simulation setting
n = 500
p = 200
rank = 10
sigma = 0.1
ratio = 0.4
rank = 10
# repetition
nrep = 20

# set your work directory to '~/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula'
source('fun_auxi.R')
source('fun_LRMC_ct.R')
source('fun_PPCA_ct.R')

info_ct = array(0, dim = c(nrep,4,2,3), 
                dimnames = list(NULL, c('LRGC','LRMC','PPCA','MIPCA'), c('LR', 'HR'), c('rate','length','time')))
for (i in 1:nrep){
  for (j in 1:2){
    if (j==1) setting = "LR" else setting = "HR"
    Xtrue = generate_LRGC(n=n, p=c(p,0,0), rank = rank, sigma = sigma, seed = i, cont.type = setting)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    
    # LRGC
    a = Sys.time()
    est = impute_mixedgc_ppca(X_obs, rank = rank, early.stop = TRUE)
    ct = ct_impute(X_obs, est)
    b = Sys.time()
    info_ct[i,1,j,1] = mean(Xtrue[loc] >= ct$lower[loc] & Xtrue[loc] <= ct$upper[loc])
    info_ct[i,1,j,2] = mean(ct$upper[loc] - ct$lower[loc])
    info_ct[i,1,j,3] = difftime(b,a,units = 'secs')
    
    # LRMC
    a = Sys.time()
    # Initial softImpute fitting
    lam0=lambda0(X_obs)
    if (j==1) index = 6 else index = 5 
    # 'index' is selected as best tuning parameter for softimpute based on imputation accuracy
    lambda=exp(seq(from=log(lam0),to=log(lam0/100),length=9))[index]
    fiti=softImpute(X_obs,lambda=lambda,rank=199,type = 'svd')
    # debiase initial imputation and construct confidence interval
    est = ct_LRMC(fiti$u, fiti$d, fiti$v, X_obs)
    b = Sys.time()
    info_ct[i,2,j,1] = mean(Xtrue[loc] >= est$lower[loc] & Xtrue[loc] <= est$upper[loc])
    info_ct[i,2,j,2] = mean(est$upper[loc] - est$lower[loc])
    info_ct[i,2,j,3] = difftime(b,a,units = 'secs')
    
    # PPCA 
    a = Sys.time()
    est = ppca_edited(X_obs, nPcs = rank) # edited only for access to noise variance estimate
    fit = list(W = est$res@loadings, sigma = est$sigma, Zimp = est$res@completeObs)
    ct = ct_impute_ppca(X_obs, fit)
    b = Sys.time()
    info_ct[i,3,j,1] = mean(Xtrue[loc] >= ct$lower[loc] & Xtrue[loc] <= ct$upper[loc])
    info_ct[i,3,j,2] = mean(ct$upper[loc] - ct$lower[loc])
    info_ct[i,3,j,3] = difftime(b,a,units = 'secs')
    
    # MIPCA
    a = Sys.time()
    est = MIPCA(X_obs, ncp = rank, nboot=100) # notice the number of samples, nboot, highly influences the performance and runtime, we follow the default setiing.
    b = Sys.time()
    imp_MI = sapply(est$res.MI, function(x){as.matrix(x)[loc]})
    ct_MI = apply(imp_MI, 1, quantile, probs = c(0.025,0.975)) # construct empirical confidence interval
    info_ct[i,4,j,1] = mean(Xtrue[loc] >= ct_MI[1,] & Xtrue[loc] <= ct_MI[2,])
    info_ct[i,4,j,2] = mean(ct_MI[2,] - ct_MI[1,])
    info_ct[i,4,j,3] = difftime(b,a,units = 'secs')
  }

  print(paste('finish iteration: ',i))
}


apply(info_ct,c(2,3,4),mean)


