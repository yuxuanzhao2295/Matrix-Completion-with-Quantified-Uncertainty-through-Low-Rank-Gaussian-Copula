# simulation about uncertainty
library(missMDA)
n = 500
p = 200
rank = 10
sigma = 0.1
ratio = 0.4
rank = 10

err_quantile_MIPCA = array(0, dim = c(20,100,2))
for (i in 1:20){
  for (j in 1:2){
    if (j==1) setting = "LR" else setting = "HR"
    Xtrue = generate_LRGC(n=n, p=c(p,0,0), rank = rank, sigma = sigma, seed = i, cont.type = setting)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    
    # 
    est = MIPCA(X_obs, ncp = 10, nboot = 20)
    imp = sapply(est$res.MI, function(x){as.matrix(x)[loc]})
    imp_var = apply(imp, 1, var)
    imp_mean = apply(imp, 1, mean)
    err_quantile_MIPCA[i,,j] = nrmse_by_prob(prob = -imp_var, xtrue = Xtrue[loc], ximp = imp_mean)
  }
  
  print(paste('finish iteration: ',i))
}
res_imp_MIPCA = list()
res_imp_MIPCA$continuous = list(err_uncertainty = err_quantile_MIPCA)

rank = 5
sigma_seq = c(0.1, 0.5)
ratio = 0.6
rank_seq = c(4,5,6)
mae_MIPCA = array(0, dim = c(20,3,2))
for (i in 1:20){
  for (j in 1:2){
    Xtrue = generate_LRGC(n=n, p=c(0,p,0), rank = rank, sigma = sigma_seq[j], seed = i, ordinal.new = TRUE)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    
    for (l in 1:3){
      est = MIPCA(X_obs, ncp = rank_seq[l], nboot =20)
      imp = sapply(est$res.MI, function(x){as.matrix(x)[loc]})
      imp = trunc.rating(apply(imp, 1, mean))
      val = Xtrue[loc]
      mae_MIPCA[i,l,j] = mean(abs(imp-val))
    }
  }
  
  print(paste('finish iteration: ',i))
}
apply(mae_MIPCA, c(2,3), mean)
apply(mae_MIPCA, c(2,3), sd)

mae_quantile_MIPCA = array(0, dim = c(20,100,2))
for (i in 1:20){
  for (j in 1:2){
    Xtrue = generate_LRGC(n=n, p=c(0,p,0), rank = rank, sigma = sigma_seq[j], seed = i, ordinal.new = TRUE)
    set.seed(i)
    loc = sample(1:prod(n*p), size = floor(prod(n*p)*ratio))
    X_obs = Xtrue
    X_obs[loc] = NA
    
    # 
    est = MIPCA(X_obs, ncp = 5, nboot = 20)
    imp = sapply(est$res.MI, function(x){as.matrix(x)[loc]})
    imp_var = apply(imp, 1, var)
    imp_mean = apply(imp, 1, mean)
    imp_mean = trunc.rating(imp_mean)
    mae_quantile_MIPCA[i,,j] = mae_by_prob(prob = -imp_var, xtrue = Xtrue[loc], ximp = imp_mean)
  }
  
  print(paste('finish iteration: ',i))
}
err = apply(mae_quantile_MIPCA, c(2,3), mean)
plot(err[,1])
plot(err[,2])

res_MIPCA_rebuttal = list(continuous = err_quantile_cont_add, ordinal = err_quantile_ordinal_add)
save(res_MIPCA_rebuttal, file = 'Rebuttal_MIPCA.RData')

err_quantile_ordinal_add = array(0, dim=c(20,100,2,2),dimnames = list(NULL, NULL, c('HighSNR',"'LowSNR"), c('LRGC','MIPCA')))
err_quantile_ordinal_add[,,,2] = mae_quantile_MIPCA
err_quantile_ordinal_add[,,,1] = err_quantile_combined_sim$ordinal[,,c(2,1),1]

p1 = plot_error_quantile_onefactor(err = err_quantile_ordinal_add[,,1,], Method = c('LRGC', 'MIPCA'),
                                   title = '1-5 Ordinal HighSNR',
                                   xlab = '',
                                   ylab = '', errorbar = TRUE)
p2 = plot_error_quantile_onefactor(err = err_quantile_ordinal_add[,,2,], Method = c('LRGC', 'MIPCA'),
                                   title = '1-5 Ordinal LowSNR',
                                   xlab = '',
                                   ylab = '', errorbar = TRUE)

err_quantile_cont_add = array(0, dim=c(20,100,2,2),dimnames = list(NULL, NULL, c('LR',"'HR"), c('LRGC','MIPCA')))
err_quantile_cont_add[,,,2] = err_quantile_MIPCA
err_quantile_cont_add[,,,1] = err_quantile_combined_sim$continuous[,,,1]
p3 = plot_error_quantile_onefactor(err = err_quantile_cont_add[,,1,], Method = c('LRGC', 'MIPCA'),
                                   title = 'Continuous LowRank',
                                   xlab = '',
                                   ylab = '', errorbar = TRUE)
p4 = plot_error_quantile_onefactor(err = err_quantile_cont_add[,,2,], Method = c('LRGC', 'MIPCA'),
                                   title = 'Continuous HighRank',
                                   xlab = '',
                                   ylab = '', errorbar = TRUE)
p3 = p3 + coord_cartesian(ylim=c(0, .8)) +  scale_color_manual(values = c("brown1", "deepskyblue"))+scale_shape_manual(values=c(15,22))
p4 = p4 + coord_cartesian(ylim=c(0, .8)) +  scale_color_manual(values = c("brown1", "deepskyblue"))+scale_shape_manual(values=c(15,22))
p1 = p1 + coord_cartesian(ylim=c(0, 1)) +  scale_color_manual(values = c("brown1", "deepskyblue"))+scale_shape_manual(values=c(15,22))
p2 = p2 + coord_cartesian(ylim=c(0, 1)) +  scale_color_manual(values = c("brown1", "deepskyblue"))+scale_shape_manual(values=c(15,22))
mylegend = g_legend(p1 +  theme(legend.position="bottom", legend.text = element_text(size=12,face="bold")))
grid.arrange(arrangeGrob(p3 + theme(legend.position = 'none'),
                         p4 + theme(legend.position = 'none'),
                         p1 + theme(legend.position = 'none'),
                         p2 + theme(legend.position = 'none'),
                         nrow=1),
             mylegend, nrow=2,heights=c(10, 1))
