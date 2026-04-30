plot.tmle <- function(plot_title,ycoord.bias=NULL,ycoord.var=NULL, ylim.bias=NULL, ylim.var=NULL,z, samplesize.label=T){

  avg.bias_ate <- get(paste0("avg.bias_ate.",z))
  bias = avg.bias_ate[avg.bias_ate$n %in% n.vec,]

  avg.variance_ate <- get(paste0("avg.variance_ate.",z))
  var = avg.variance_ate[avg.variance_ate$n %in% n.vec,]

  true_asyvar = VAR.ATE[z]

  yaxis_bias = TeX(r'($\sqrt{n}$-Bias)')
  yaxis_var = TeX(r'($n$-Variance)')


  ######################
  # bias plot
  ######################

  p.bias <- ggplot(data=bias, aes(x=n, y=bias)) +
    geom_line(linewidth=0.7)+
    geom_hline(yintercept=0, linetype="dashed",color="blue",linewidth=0.8)+
    coord_cartesian(ylim=ycoord.bias)+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=13,face = "bold")) +
    theme(axis.text.x = element_text(angle = -45))+
    scale_x_continuous(breaks=n.vec)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=25,
                  position=position_dodge(0.05))+
    geom_point(size=3,color="red")+
    labs(title="",x="", y = yaxis_bias)

  if (!is.null(ylim.bias)){

    p.bias <- p.bias + ylim(ylim.bias[1], ylim.bias[2])

  }

  ######################
  # variance plot
  ######################

  p.variance <- ggplot(data=var, aes(x=n, y=variance)) +
    geom_line(linewidth=0.7)+
    geom_hline(yintercept=true_asyvar, linetype="dashed",color="blue",linewidth=0.8)+
    coord_cartesian(ylim=ycoord.var)+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=13,face="bold")) +
    theme(axis.text.x = element_text(angle = -45))+
    scale_x_continuous(breaks=n.vec)+
    geom_errorbar(aes(ymin=lower, ymax=upper), width=25,
                  position=position_dodge(0.05),linewidth=0.5)+
    geom_point(size=3,color="red")+
    labs(title="",x=ifelse(samplesize.label,'Sample size n',''), y = yaxis_var)

  if (!is.null(ylim.var)){
   p.variance <- p.variance + ylim(ylim.var[1], ylim.var[2])
  }

  ######################
  # arrangement
  ######################

  # arrange multiple plots
  p <- plot_grid(
    p.bias
    ,p.variance
    , align = "hv"
    , ncol = 1
  )
  p <- annotate_figure(p, left = text_grob(TeX(plot_title), color = "black", face = "bold", size = 17,rot = 90))

  return(p)
}
