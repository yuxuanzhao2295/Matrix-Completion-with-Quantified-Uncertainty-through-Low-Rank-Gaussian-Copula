# this script was ran on cluster in parallelization with 10 cores due to long time searching
library(doParallel)
registerDoParallel(20)

library(mixedgcImp)
library(softImpute)


load("~/impute/movielens_50obs/moviesWith50ratings.RData")

ratiosoftImpute = 1/4
res = foreach(i = 1:5)%dopar%{
  split = data_split(movie.50obs, seed = i)
  loc = split$loc_test
  loc = sort(loc)
  info = array(0, dim = c(100, 2))
  
  xc = biScale(split$Xtrain,row.scale = FALSE,col.scale = FALSE, trace = FALSE)
  lam0=lambda0(xc)  
  fiti = softImpute(xc,lambda=lam0*ratiosoftImpute, rank.max = 500, type = 'svd')
  ximp = trunc.rating(complete(xc, fiti))
  
  r = varest_softImpute(X = split$Xtrain,  loc = loc, ratio = ratiosoftImpute, parallel = TRUE, maxrank = 500, seed = i)
  # smaller variance, larger reliabilisty
  info[,1] = mae_by_reliability(r = -r, xtrue = movie.50obs[loc], ximp = ximp[loc])
  info[,2] = rmse_by_reliability(r = -r, xtrue = movie.50obs[loc], ximp = ximp[loc])
  
  info
}
resReliability_softimpute_movielens = array(0, dim = c(5,100,2), 
                                            dimnames = list(NULL,NULL,c('MAE','RMSE')))
for (i in 1:5) resReliability_softimpute_movielens[i,,] = res[[i]]
save(resReliability_softimpute_movielens, file = 'resReliability_softimpute_movielens.RData')
