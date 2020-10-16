# simulation in producing confidence interval
library(mixedgcImp)
library(softImpute)
library(pcaMethods)
# simulation setting
n = 500
p = 200
rank = 5
sigmaSeq = c(0.1,0.5)
ratio = 0.6
nrep = 20 # repetition

# search sequence
rankLRGC = 3:11
rankPPCA = 3:11
ratiosoftImpute = exp(seq(from=log(1),to=log(1/100),length=9))


# set your work directory to '~/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula'
source('fun_auxi.R')

resTuning_ord = array(0, dim = c(nrep,9,2,2,3), 
                       dimnames = list(NULL, NULL, c('MAE','time'), c('HighSNR', 'LowSNR'),c('LRGC','PPCA','softImpute')))
for (i in 1:nrep){
  for (j in 1:2){
    Xtrue = generate_LRGC(n=n, p=c(0,p,0), rank = rank, sigma = sigmaSeq[j], seed = i)
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
      resTuning_ord[i,s,1,j,1] = cal_mae(est$Ximp, X_obs, Xtrue)
      resTuning_ord[i,s,2,j,1] = difftime(b,a,units = 'secs')
      
      # PPCA
      a = Sys.time()
      est = pca(X_obs, method = 'ppca', nPcs = rankPPCA[s])
      b = Sys.time()
      resTuning_ord[i,s,1,j,2] = cal_mae(trunc.rating(est@completeObs), X_obs, Xtrue)
      resTuning_ord[i,s,2,j,2] = difftime(b,a,units = 'secs')
      
      # softImpute
      a = Sys.time()
      fiti=softImpute(xc,lambda=lam0 * ratiosoftImpute[s],rank=199, warm=warm, type = 'svd')
      ximp = complete(xc, fiti, unscale = TRUE)
      b = Sys.time()
      resTuning_ord[i,s,1,j,3] = cal_mae(trunc.rating(ximp), X_obs, Xtrue)
      resTuning_ord[i,s,2,j,3] = difftime(b,a,units = 'secs')
      warm=fiti
    }
  }
  print(paste('finish iteration: ',i))
}

save(resTuning_ord, file = 'ResTuning_LRGC_PPCA_softImpute_SimOrdinal.RData')
err = apply(resTuning_ord, c(2,3,4,5), mean)
