library(mixedgcImp)
library(softImpute)
library(pcaMethods)
# simulation setting
n = 500
p = 200
rank = 10
sigma = 0.1
ratio = 0.4
nrep = 20 # repetition

# search sequence
rankLRGC = c(6,7,8,9,10,11,12,13,14)
rankPPCA = c(6,7,8,9,10,11,12,13,14)
ratiosoftImpute = exp(seq(from=log(1),to=log(1/100),length=9))


# set your work directory to '~/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula'
source('fun_auxi.R')

resTuning_cont = array(0, dim = c(nrep,9,2,2,3), 
            dimnames = list(NULL, NULL, c('NRMSE','time'), c('LR', 'HR'),c('LRGC','PPCA','softImpute')))
for (i in 1:nrep){
  for (j in 1:2){
    if (j==1) setting = "LR" else setting = "HR"
    Xtrue = generate_LRGC(n=n, p=c(p,0,0), rank = rank, sigma = sigma, seed = i, cont.type = setting)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    
    #preparation for softimpute
    xc = biScale(X_obs,row.scale=FALSE,col.scale=FALSE)
    lam0=lambda0(xc)
    warm=NULL
    
    # start search over the grid
    for (s in 1:9){
      # LRGC
      a = Sys.time()
      est = impute_mixedgc_ppca(X_obs, rank = rankLRGC[s]) 
      b = Sys.time()
      resTuning_cont[i,s,1,j,1] = cal_rmse(est$Ximp, X_obs, Xtrue)
      resTuning_cont[i,s,2,j,1] = difftime(b,a,units = 'secs')
      
      # PPCA
      a = Sys.time()
      est = pca(X_obs, method = 'ppca', nPcs = rankPPCA[s])
      b = Sys.time()
      resTuning_cont[i,s,1,j,2] = cal_rmse(est@completeObs, X_obs, Xtrue)
      resTuning_cont[i,s,2,j,2] = difftime(b,a,units = 'secs')
      
      # softImpute
      a = Sys.time()
      fiti=softImpute(xc,lambda=lam0 * ratiosoftImpute[s],rank=199, warm=warm, type = 'svd')
      ximp = complete(xc, fiti, unscale = TRUE)
      b = Sys.time()
      resTuning_cont[i,s,1,j,3] = cal_rmse(ximp, X_obs, Xtrue)
      resTuning_cont[i,s,2,j,3] = difftime(b,a,units = 'secs')
      warm=fiti
    }
  }
  print(paste('finish iteration: ',i))
}
save(resTuning_cont, file = 'ResTuning_LRGC_PPCA_softImpute_SimContinuous.RData')
err = apply(resTuning_cont, c(2,3,4,5), mean)
