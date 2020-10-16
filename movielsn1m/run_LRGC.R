setwd('/Users/yuxuan/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula')
load("movielsn1m/moviesWith50ratings.RData")
library(mixedgcImp)
nrep = 5
nquantiles = 100

rankCopula = 10
errAll = array(0, dim = c(nrep,4), dimnames = list(NULL, c('MAE','RMSE','RunTime','ReliabiliyRuntime')))
errReliability = array(0, dim = c(nrep,nquantiles,2), dimnames = list(NULL, NULL, c('MAE','RMSE')))
for (i in 1:nrep){
  split = data_split(movie.50obs, seed = i)
  
  # 
  a = Sys.time()
  est = impute_mixedgc_ppca(split$Xtrain, rank = rankCopula, eps = 1e-4)
  b = Sys.time()
  prob = prob_impute(split$Xtrain, est)$prob
  c = Sys.time()
  
  e_test = est$Ximp[split$loc_test] - movie.50obs[split$loc_test] # test error
  errAll[i,1] = mean(abs(e_test))
  errAll[i,2] = sqrt(mean(e_test^2))
  errAll[i,3] = difftime(b,a,units = 'mins')
  errAll[i,4] = difftime(c,b,units = 'mins')
  #
  prob = prob_impute(split$Xtrain, est)$prob
  loc_test = split$loc_test
  errReliability[i,,1] = mae_by_prob(prob[loc_test],movie.50obs[loc_test],est$Ximp[loc_test], l=nquantiles)
  errReliability[i,,2] = rmse_by_prob(prob[loc_test],movie.50obs[loc_test],est$Ximp[loc_test], l=nquantiles)
  
  print(paste('finish repetition ',i))
}

apply(errAll,2,mean)
e = apply(errReliability, c(2,3),mean)
plot(rev(e[,1])) 
plot(rev(e[,2]))

# save results if needed
ResLRGC_movielens = list(errAll = errAll, errReliability = errReliability)
save(ResLRGC_movielens, file='ResLRGC_movielens.RData')
