setwd('~/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula')
load("GLRM_julia/Res_sim_GLRM.RData")
load("simulation/ResTuning_LRGC_PPCA_softImpute_SimBinary.RData")
load("simulation/ResTuning_LRGC_PPCA_softImpute_SimOrdinal.RData")
load("simulation/ResTuning_LRGC_PPCA_softImpute_SimContinuous.RData")
library(gridExtra)

err_allsim_cont = array(0, dim=c(20,9,2,5),
                                   dimnames = list(NULL, Tuning = NULL, Setting = c('LowRank','HighRank'), Method = c('LRGC','PPCA','softImpute','GLRM-l2','MMC')))
err_allsim_ord = array(0, dim=c(20,9,2,5),
                                  dimnames = list(NULL, Tuning = NULL, Setting = c('HighSNR','LowSNR'), Method = c('LRGC','PPCA','softImpute','GLRM-BvS', 'GLRM-l1')))
err_allsim_bin = array(0, dim=c(20,9,2,5),
                                  dimnames = list(NULL, Tuning = NULL, Setting = c('HighSNR','LowSNR'), Method = c('LRGC','PPCA','softImpute', 'GLRM-logistic','GLRM-hinge')))
err_allsim_cont[,,,1:3] = resTuning_cont[,,1,,]
err_allsim_cont[,,,4] = Res_sim_MMMF$resTuning_cont[,,1,]
# load from old data, to be updated
load("/Users/yuxuan/Box/imputation/simulation_low_rank/Continuous_MMC.RData")
err_allsim_cont[,,1,5] = info_MMC_LR[,,1,1]
err_allsim_cont[,,2,5] = info_MMC_HR[,,1,1]

err_allsim_ord[,,,1:3] = resTuning_ord[,,1,,]
err_allsim_ord[,,,4:5] = Res_sim_MMMF$resTuning_ord[,,1,,]

err_allsim_bin[,,,1:3] = resTuning_bin[,,1,,]
err_allsim_bin[,,,4:5] = Res_sim_MMMF$resTuning_bin[,,1,,]

# continuous 
titles = dimnames(err_allsim_cont)[[4]]
xnames = c('rank', 'rank', 'penalization', 'penalization', 'stepsize')
xticks = list(6:14, 6:14, 
             round(seq(from=log(1),to=log(1/100),length=9),1),
             round(seq(from=log(1/4),to=log(1/100),length=9),1),
             seq(3,19,2)
             )
p.contlist = list()
for (i in 1:5){
  if (i == 1) ylab = "NRMSE" else yalb = ""
  p.contlist[[i]] = plot_imp(err_allsim_cont[,,,i], xtick = xticks[[i]], xlab=xnames[i], ylab=ylab, title = titles[i], ylims = c(0,1))
}

# ordinal
titles = dimnames(err_allsim_ord)[[4]]
xnames = c('rank', 'rank', 'penalization', 'penalization', 'penalization')
xticks = list(3:11, 3:11, 
              round(seq(from=log(1),to=log(1/100),length=9),1),
              round(seq(from=log(1/4),to=log(1/100),length=9),1),
              round(seq(from=log(1/4),to=log(1/100),length=9),1))
p.ordlist = list()
for (i in 1:5){
  if (i == 1) ylab = "MAE" else yalb = ""
  p.ordlist[[i]] = plot_imp(err_allsim_ord[,,,i], xtick = xticks[[i]], xlab=xnames[i], ylab=ylab, title = titles[i], ylims = c(0,1))
}

# binary
titles = dimnames(err_allsim_bin)[[4]]
xnames = c('rank', 'rank', 'penalization', 'penalization', 'penalization')
xticks = list(3:11, 3:11, 
              round(seq(from=log(1),to=log(1/100),length=9),1),
              round(seq(from=log(1),to=log(1/100),length=9),1),
              round(seq(from=log(1),to=log(1/100),length=9),1))
p.binlist = list()
for (i in 1:5){
  if (i == 1) ylab = "MAE" else yalb = ""
  p.binlist[[i]] = plot_imp(err_allsim_bin[,,,i], xtick = xticks[[i]], xlab=xnames[i], ylab=ylab, title = titles[i], ylims = c(0,0.3))
}

mylegend1 = g_legend(p.contlist[[1]] +  theme(legend.position="bottom"))
mylegend2 = g_legend(p.ordlist[[1]] +  theme(legend.position="bottom"))

grid.arrange(arrangeGrob(p.contlist[[1]],p.contlist[[2]],p.contlist[[3]],p.contlist[[4]],p.contlist[[5]],
                         nrow=1),
             mylegend1,
             arrangeGrob(p.ordlist[[1]],p.ordlist[[2]],p.ordlist[[3]],p.ordlist[[4]],p.ordlist[[5]],
                         p.binlist[[1]],p.binlist[[2]],p.binlist[[3]],p.binlist[[4]],p.binlist[[5]],
                         nrow=2),
             mylegend2, 
             nrow=4,heights=c(10, 1,20,1))
