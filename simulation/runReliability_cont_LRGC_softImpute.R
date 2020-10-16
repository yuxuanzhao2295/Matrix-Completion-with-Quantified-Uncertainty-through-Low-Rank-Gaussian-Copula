# simulation in producing confidence interval
library(mixedgcImp)
library(softImpute)
source('fun_auxi.R')
# simulation setting
n = 500
p = 200
rank = 10
sigma = 0.1
ratio = 0.4
nrep = 20 # repetition
nquantiles = 100

# optimal tuning parameter under two settings LR, HR
rankLRGC = c(10, 10)
ratiosoftImpute = exp(seq(from=log(1),to=log(1/100),length=9))[c(6,5)]
nfolds = 10


errReliability_cont = array(0, dim = c(nrep,nquantiles,2,2))
for (i in 1:20){
  for (j in 1:2){
    if (j==1) setting = "LR" else setting = "HR"
    Xtrue = generate_LRGC(n=n, p=c(p,0,0), rank = rank, sigma = sigma, seed = i, cont.type = setting)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    loc = which(is.na(X_obs))
    
    # LRGC
    est = impute_mixedgc_ppca(X_obs, rank = rankLRGC[j]) 
    r = reliability_cont(X_obs, est)
    errReliability_cont[i,,j,1] = nrmse_by_reliability(r = r[loc], xtrue = Xtrue[loc], ximp = est$Ximp[loc])
    
    # softImpute
    xc = biScale(X_obs,row.scale=FALSE,col.scale=FALSE) # initial fit
    lam0=lambda0(xc)*ratiosoftImpute[j]
    fiti=softImpute(xc,lambda=lam0,rank=199,type = 'svd')
    ximp = complete(xc, fiti)
    r = varest_softImpute(X_obs, ratio = ratiosoftImpute[j], seed = i)
    # smaller variance, larger reliabilisty
    errReliability_cont[i,,j,2] = nrmse_by_reliability(r = -r, xtrue = Xtrue[loc], ximp = ximp[loc])
  }

  print(paste('finish iteration: ',i))
}
save(errReliability_cont, file = 'ResReliability_LRGC_softImpute_SimContinuous.RData')

err = apply(errReliability_cont, c(2,3,4), mean)
plot(err[,1,1], col ='red', ylim = c(0,1))
points(err[,1,2], col ='green', ylim = c(0,1))
plot(err[,2,1], col ='red', ylim = c(0,1))
points(err[,2,2], col ='green', ylim = c(0,1))
