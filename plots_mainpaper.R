library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)
setwd('~/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula')
source('fun_plots.R')

# simulation plot Figure 1

errReliability_allsim_cont = array(0, dim=c(20,100,2,4),
                                dimnames = list(NULL, Quantile = 1:100, Setting = c('LowRank','HighRank'), Method = c('LRGC','MI+softImpute','MI+GLRM-l2','MI-PCA')))
errReliability_allsim_ord = array(0, dim=c(20,100,2,4),
                               dimnames = list(NULL, Quantile = 1:100, Setting = c('HighSNR','LowSNR'), Method = c('LRGC','MI+softImpute','MI+GLRM-BvS', 'MI-PCA')))
errReliability_allsim_bin = array(0, dim=c(20,100,2,4),
                               dimnames = list(NULL, Quantile = 1:100, Setting = c('HighSNR','LowSNR'), Method = c('LRGC','MI+softImpute', 'MI+GLRM-logistic', 'MI+GLRM-hinge')))

setwd("/Users/yuxuan/Documents/GitHub/Matrix-Completion-with-Quantified-Uncertainty-through-Low-Rank-Gaussian-Copula")
load("simulation/ResReliability_LRGC_softImpute_SimBinary.RData")
load("simulation/ResReliability_LRGC_softImpute_SimContinuous.RData")
load("simulation/ResReliability_LRGC_softImpute_SimOrdinal.RData")
load("simulation/Rebuttal_MIPCA.RData")
load("GLRM_julia/Res_sim_GLRM.RData")
load("GLRM_julia/errReliability_bin_GLRM_hinge.RData")

errReliability_allsim_cont[,,,1:2] = errReliability_cont
errReliability_allsim_cont[,,,3] = Res_sim_MMMF$errReliability_cont_MMMF
errReliability_allsim_cont[,,,4] = res_MIPCA_rebuttal$continuous[,,,2]

errReliability_allsim_ord[,,,1:2] = errReliability_ord
errReliability_allsim_ord[,,,3] = Res_sim_MMMF$errReliability_ord_MMMF
errReliability_allsim_ord[,,,4] = res_MIPCA_rebuttal$ord[,,,2]

errReliability_allsim_bin[,,,1:2] = errReliability_bin
errReliability_allsim_bin[,,,3] = Res_sim_MMMF$errReliability_bin_MMMF
errReliability_allsim_bin[,,,4] = errReliability_bin_MMMF_hinge 


p1 = plot_error_quantile_twofactor(err = errReliability_allsim_cont, title = 'Continuous',xlab = 'Percentage of entries selected',ylab = 'NRMSE', errorbar = TRUE,
                                   ylims=c(0,0.8), colorvals = c("brown1", "deepskyblue", "deeppink",'darkolivegreen3'), shapevals = c(15,18,3,4))
p2 = plot_error_quantile_twofactor(err = errReliability_allsim_ord, title = '1-5 Ordinal',xlab = 'Percentage of entries selected',ylab = 'MAE', errorbar = TRUE,
                                   ylims=c(0,1), colorvals = c("brown1", "deepskyblue", "darkorange",'darkolivegreen3'), shapevals = c(15,18,17,4))
p3 = plot_error_quantile_twofactor(err = errReliability_allsim_bin, title = 'Binary',xlab = 'Percentage of entries selected',ylab = 'MAE', errorbar = TRUE,
                                   ylims=c(0,0.25), colorvals = c("brown1", "deepskyblue", "cornsilk4","mediumpurple"), shapevals = c(15,18,16,11))
mylegend = produce_legend(names = c("LRGC", "MI+softImpute", "MI+GLRM-l2", "MI-PCA", "MI+GLRM-BvS", "MI+GLRM-logistic", "MI+GLRM-hinge"),
                          colorvals =  c("brown1", "deepskyblue", "deeppink",'darkolivegreen3', 'darkorange', 'cornsilk4', 'mediumpurple'), 
                          shapevals = c(15,18,3,4,17,16,11))

grid.arrange(arrangeGrob(p1,p2,p3,nrow=1),
             mylegend, nrow=2,heights=c(9, 1))


# movielens result figure 2
errReliability_allmovielens = array(0, dim=c(5,100,2,3),
                                   dimnames = list(NULL, Quantile = 1:100, Metric = c('MAE','RMSE'), Method = c('LRGC','MI+softImpute','MI+GLRM-l2')))

load("GLRM_julia/Res_movielens_GLRM.RData")
load("movielsn1m/resReliability_softimpute_movielens.RData")
#load("movielsn1m/ResLRGC_movielens.RData") # corrupted data
errReliability_allmovielens[,,,1] = ResLRGC_movielens$errReliability
errReliability_allmovielens[,,,2] = resReliability_softimpute_movielens
errReliability_allmovielens[,,,3] = res_MMMF$errReliability

p4 = plot_error_quantile_onefactor(err = errReliability_allmovielens[,,1,], title = '',xlab = '',ylab = 'MAE', errorbar = TRUE,
                                   ylims=c(0,1), colorvals = c("brown1", "deepskyblue", "darkorange"), shapevals = c(15,18,17))

# plot from history data, to be updated 
load("/Users/yuxuan/Downloads/Uncertainty_movie50obs.RData")
load("/Users/yuxuan/Jupyter_notebook/Uncertainty_movie50obs_MMMF.RData")
err = array(0,dim = c(5,100,3),dimnames = list(NULL,Quantile=1:100,Method=c('Our: LRGC','MI+softImpute','MI+GLRM-BvS')))
err[,,1] = error_quantile_movie50$EM[,,1]
err[,,2] = error_quantile_movie50$soft[,,1]
err[,,3] = info_q[,,1]
p4 = plot_error_quantile_onefactor(err = err, title = 'MovieLens 1M',xlab = 'Percentage of Entries Selected',ylab = 'MAE', errorbar = TRUE,
                                   ylims=c(0,0.75), colorvals = c("brown1", "deepskyblue", "darkorange"), shapevals = c(15,18,17), position="right")

err = array(0,dim = c(5,100,3),dimnames = list(NULL,Quantile=1:100,Method=c('LRGC','10fold-softImpute','10fold-MMMF-BvS')))
err[,,1] = sqrt(error_quantile_movie50$EM[,,2])
err[,,2] = sqrt(error_quantile_movie50$soft[,,2])
err[,,3] = info_q[,,2]
p5 = plot_error_quantile_onefactor(err = err, title = '',xlab = '',ylab = 'RMSE', errorbar = TRUE,
                                   ylims=c(0,1), colorvals = c("brown1", "deepskyblue", "darkorange"), shapevals = c(15,18,17))
mylegend = produce_legend(names = c("LRGC", "MI+softImpute", "MI+GLRM-BvS"),
                          colorvals =  c("brown1", "deepskyblue", 'darkorange'), 
                          shapevals = c(15,18,17), position = "right")
grid.arrange(arrangeGrob(p4,p5,nrow=1),
             mylegend, ncol=2,widths=c(3, 1),
             top = textGrob("Movielens 1M",gp=gpar(fontsize=10, font=2)),
             bottom = textGrob("Percentages of entries selected",gp=gpar(fontsize=10, font=2)))
