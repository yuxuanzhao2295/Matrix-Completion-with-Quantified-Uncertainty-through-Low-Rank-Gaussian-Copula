# this script was ran on cluster in parallelization with 10 cores due to long time searching
library(doParallel)
registerDoParallel(10)
getDoParWorkers()

library(mixedgcImp)
library(softImpute)

setwd('/Users/yuxuan/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula')
load("movielsn1m/moviesWith50ratings.RData")

rankLRGC  = c(8, 10, 12, 14)
ratiosoftImpute = c(1/2, 1/4, 1/6, 1/8)
res = foreach(i = 1:10)%dopar%{
  if (i<=5){
    # first 5 experiments for copula_EM
    split = data_split(movie.50obs, seed = i)
    
    info = array(0, dim = c(4,4))
    for (j in 1:4){
      est = impute_mixedgc_ppca(split$Xtrain, rank = rankLRGC[j], eps = 1e-4)
      e_test = est$Ximp[split$loc_test] - movie.50obs[split$loc_test]
      e_val = est$Ximp[split$loc_val] - movie.50obs[split$loc_val]
      info[j,1] = mean(abs(e_test))
      info[j,2] = sqrt(mean(e_test^2))
      info[j,3] = mean(abs(e_val))
      info[j,4] = sqrt(mean(e_val^2))
    }
    info
  }else{
    # next 5 experiments for softimpute
    split = data_split(movie.50obs, seed = i-5)
    
    xc = biScale(split$Xtrain,row.scale = FALSE,col.scale = FALSE, trace = FALSE)
    lam0=lambda0(xc) 
    info = array(0, dim = c(4,6))
    warm = NULL
    for (j in 1:4){
      fiti = softImpute(xc,lambda=lam0*ratiosoftImpute[j], rank.max = 500, type = 'svd', warm = warm)
      xhat = trunc.rating(complete(xc, fiti))
      warm = fiti
      e_test = xhat[split$loc_test] - movie.50obs[split$loc_test]
      e_val = xhat[split$loc_val] - movie.50obs[split$loc_val]
      info[j,1] = mean(abs(e_test))
      info[j,2] = sqrt(mean(e_test^2))
      info[j,3] = mean(abs(e_val))
      info[j,4] = sqrt(mean(e_val^2))
    }
    info
  }
}

resTuning_movielens = array(0, dim = c(5,4,4,2), 
                            dimnames = list(NULL, NULL,
                                            c('MSE: test', 'RMSE: test','MSE: validation', 'RMSE:validation'), 
                                            c('LRGC','softImpute')))

for (i in 1:5){
  resTuning_movielens[i,,1] = res[[i]]
  resTuning_movielens[i,,2] = res[[i+5]]
}