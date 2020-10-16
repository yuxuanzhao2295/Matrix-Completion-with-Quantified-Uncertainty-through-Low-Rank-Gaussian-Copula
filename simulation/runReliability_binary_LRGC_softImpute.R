# simulation in producing confidence interval
library(mixedgcImp)
library(softImpute)
# set your work directory to '~/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula'
source('fun_auxi.R')
# simulation setting
n = 500
p = 200
rank = 5
sigmaSeq = c(0.1,0.5)
ratio = 0.6
nrep = 20 # repetition
nquantiles = 100

# search sequence
rankLRGC = c(5,5)
ratiosoftImpute = exp(seq(from=log(1),to=log(1/100),length=9))[c(4,3)]

errReliability_bin = array(0, dim = c(nrep,nquantiles,2,2))
for (i in 1:nrep){
  for (j in 1:2){
    Xtrue = generate_LRGC(n=n, p=c(0,0,p), rank = rank, sigma = sigmaSeq[j], seed = i)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    loc = which(is.na(X_obs))
    
    # LRGC
    est = impute_mixedgc_ppca(X_obs, rank = rankLRGC[j]) 
    r = reliability_ord(X_obs, est)
    errReliability_bin[i,,j,1] = mae_by_reliability(r = r[loc], xtrue = Xtrue[loc], ximp = est$Ximp[loc])
    
    # softImpute
    xc = biScale(X_obs,row.scale=FALSE,col.scale=FALSE) # initial fit
    lam0=lambda0(xc)*ratiosoftImpute[j]
    fiti=softImpute(xc,lambda=lam0,rank=199,type = 'svd')
    ximp = trunc.rating(complete(xc, fiti), xmin = 0, xmax = 1)
    r = varest_softImpute(X_obs, ratio = ratiosoftImpute[j], seed = i)
    # smaller variance, larger reliabilisty
    errReliability_bin[i,,j,2] = mae_by_reliability(r = -r, xtrue = Xtrue[loc], ximp = ximp[loc])
  }
  print(paste('finish iteration: ',i))
}
save(errReliability_bin, file = 'ResReliability_LRGC_softImpute_SimBinary.RData')


err = apply(errReliability_bin, c(2,3,4), mean)
plot(err[,1,1], col ='red', ylim = c(0,0.25))
points(err[,1,2], col ='green')
plot(err[,2,1], col ='red', ylim = c(0,0.25))
points(err[,2,2], col ='green')
