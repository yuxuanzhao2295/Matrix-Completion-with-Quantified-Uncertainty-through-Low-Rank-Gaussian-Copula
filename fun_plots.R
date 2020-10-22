plot_imp = function(err, xtick, xlab="", ylab ="", title ="", colorname = NULL, position = 'none', ylims = NULL){
  # err is a 3-dimensional array
  require(ggplot2)
  require(reshape2)
  dimnames(err)[[2]] = xtick
  error = apply(err, c(2,3),mean)
  error = melt(error)
  colnames(error)[1] = xlab
  colnames(error)[3] = "mean"
  col_names = colnames(error) 
  error$sd = as.vector(apply(err, c(2,3),sd))
  error[,2] = as.factor(error[,2])
  error$ymin = error[,3] - error$sd
  error$ymax = error[,3] + error$sd
  
  
  p = ggplot(error, aes_string(x=col_names[1], y=col_names[3], group = col_names[2], colour = col_names[2], shape = col_names[2])) +
    geom_point(size =3) +
    geom_line()  +
    geom_errorbar(aes(ymin=ymin, ymax=ymax, width=.005)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))+
    theme_classic() +
    theme(axis.text=element_text(size=13, face = "bold"),
          axis.text.x = element_text(size=11,face="bold"),
          axis.text.y = element_text(size=11,face="bold"),
          axis.title=element_text(size=13,face="bold"),
          plot.title = element_text(hjust = 0.5, size=14,face="bold"),
          legend.text = element_text(size=12,face="bold"),
          legend.position = position) 
  p = p + scale_x_continuous(col_names[1],breaks = xtick,labels = as.character(xtick)) 
  if (!is.null(colorname)) p = p + scale_fill_manual(values=colorname)
  if (!is.null(ylims)) p = p + coord_cartesian(ylim=ylims)
  p + labs(x = xlab, y=ylab)
}


plot_error_quantile_twofactor = function(err, col_names, title, xlab, ylab, errorbar=TRUE,ylims=NULL,colorvals,shapevals, position = "None"){
  require(ggplot2)
  require(reshape2)
  # dim: rep * quantile * setting * methods
  index = seq(10,100,by = 10)
  error = apply(err[,index,,],c(2,3,4),mean)
  for (i in 1:dim(error)[2]){
    for (j in 1:dim(error)[3]){
      error[,i,j] = rev(error[,i,j])
    }
  }
  error = melt(error)
  colnames(error)[4] = 'mean'
  error$sd = as.vector(apply(err[,index,,], c(2,3,4),sd))
  error$ymin = error[,4] - error$sd
  error$ymax = error[,4] + error$sd
  #error$group = as.factor(paste(error[,2],error[,3]))
  #col_names = c(col_names,'group')
  col_names = colnames(error)
  
  
  p = ggplot(error, 
             aes_string(x=col_names[1], y=col_names[4], 
                        group = col_names[3], colour = col_names[3], shape = col_names[3])) +
    geom_point(size=3) +
    geom_line()  +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))+
    #facet_grid(Setting ~., scales = "free")+
    facet_grid(.~Setting, scales = "free")+
    labs(x=xlab, y = ylab)+
    theme_classic() +
    theme(axis.text=element_text(size=14, face = "bold"),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(size=14,face="bold"),
          axis.text.y = element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5, size=14,face="bold"),
          strip.text.x = element_text(size = 14, face = "bold"),
          legend.text = element_text(size=15,face="bold"),
          legend.position = position,
          plot.margin=unit(c(0.1,0.1,0.5,0.1), "cm"))
  if (errorbar) p = p + geom_errorbar(aes(ymin=ymin, ymax=ymax, width=.005))
  p = p + scale_x_continuous(breaks = c(1,25,50,75,100)) 
  if (!is.null(ylims)) p = p + coord_cartesian(ylim=ylims) 
  p + scale_color_manual(values = colorvals)+scale_shape_manual(values=shapevals)
}

plot_error_quantile_onefactor = function(err, xlab='', ylab='', title='',errorbar=TRUE, ylims=NULL,colorvals,shapevals, position = "None"){
  require(ggplot2)
  require(reshape2)
  # dim: rep * quantile * rank
  index = seq(10,100,by = 10)
  error = apply(err[,index,],c(2,3),mean)
  for (i in 1:dim(error)[2]) error[,i] = rev(error[,i])
  error = melt(error)
  colnames(error)[3] = 'mean'
  error$sd = as.vector(apply(err[,index,], c(2,3),sd))
  error$ymin = error[,3] - error$sd
  error$ymax = error[,3] + error$sd
  col_names = colnames(error)
  
  p = ggplot(error, 
             aes_string(x=col_names[1], y=col_names[3], 
                        group = col_names[2], colour = col_names[2], shape = col_names[2]))  +
    geom_point(size = 3) +
    geom_line()  +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))+
    #facet_grid(.~Setting, scales = "free")+
    labs(x=xlab, y = ylab)+
    theme_classic() +
    theme(axis.text=element_text(size=14, face = "bold"),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x = element_text(size=14,face="bold"),
          axis.text.y = element_text(size=14,face="bold"),
          plot.title = element_text(hjust = 0.5, size=14,face="bold"),
          legend.text = element_text(size=12,face="bold"),
          legend.position = position,
          plot.margin=unit(c(0.1,0.2,0,0.2), "cm"))
  if (errorbar) p = p + geom_errorbar(aes(ymin=ymin, ymax=ymax, width=.005))
  p + scale_x_continuous(breaks = c(1,25,50,75,100)) 
  if (!is.null(ylims)) p = p + coord_cartesian(ylim=ylims) 
  p + scale_color_manual(values = colorvals)+scale_shape_manual(values=shapevals)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

produce_legend <- function(names, colorvals, shapevals, position="bottom"){
  l = length(names)
  x =array(0, dim=c(20,100,l),dimnames = list(NULL, Quantile = 1:100, Method = names))
  x[,,1:l] = rnorm(prod(dim(x)))
  p = plot_error_quantile_onefactor(err=x, colorvals=colorvals, shapevals=shapevals)# + guides(colour = guide_legend(nrow = 1))
  p = p + theme(plot.margin=unit(c(1.5,0,-0.5,0), "cm"))
  g_legend(p +  theme(legend.position=position, legend.text = element_text(size=10,face="bold"), legend.key.size = unit(2,"line")))
}

