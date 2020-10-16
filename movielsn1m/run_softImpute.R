setwd('/Users/yuxuan/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula')
load("movielsn1m/moviesWith50ratings.RData")
library(softImpute)
nrep = 5
nquantiles = 100

errAll = array(0, dim = c(nrep,4,2))
errReliability = array(0, dim = c(nrep,nquantiles,2,2))
for (i in 1:nrep){
  split = data_split(movie.50obs, seed = i)
  
  # 
  a = Sys.time()
  est = impute_mixedgc_ppca(split$Xtrain, rank = rankCopula, eps = 1e-4)
  b = Sys.time()
  prob = prob_impute(split$Xtrain, est)$prob
  c = Sys.time()
  
  e_test = est$Ximp[split$loc_test] - movie.50obs[split$loc_test] # test error
  #e_val = est$Ximp[split$loc_val] - movie.50obs[split$loc_val] # validation error
  errAll[i,1,1] = mean(abs(e_test))
  errAll[i,2,1] = sqrt(mean(e_test^2))
  errAll[i,3,1] = difftime(b,a,units = 'mins')
  errAll[i,4,1] = difftime(c,b,units = 'mins')
  #
  prob = prob_impute(split$Xtrain, est)$prob
  loc_test = split$loc_test
  errReliability[i,,1,1] = mae_by_prob(prob[loc_test],movie.50obs[loc_test],est$Ximp[loc_test], l=nquantiles)
  errReliability[i,,2,1] = rmse_by_prob(prob[loc_test],movie.50obs[loc_test],est$Ximp[loc_test], l=nquantiles)
  
  print(paste('finish repetition ',i))
}

apply(errAll,c(2,3),mean)
e = apply(errReliability, c(2,3,4),mean)[,,1]
plot(rev(e[,1])) 
plot(rev(e[,2]))