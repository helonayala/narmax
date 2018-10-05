plot_xcorrel = function(e,u) {
  
  source("multiplot.R")
  source("crossco.R")
  library(reshape2)
  library(latex2exp)
  
  maxlag = 25
  N = length(u)
  
  conf_factor = 1.96/sqrt(N)
  lag_vec = -maxlag:maxlag
  
  EE = crossco(e,e,maxlag)
  UE = crossco(u,e,maxlag)
  EEU = crossco(e[1:(N-1)],e[2:N]*u[2:N],maxlag)
  U2E = crossco(u^2 - mean(u^2),e,maxlag)
  U2E2 = crossco(u^2 - mean(u^2),e^2,maxlag)
  
  EEdf = melt(data.frame(EE = EE,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  UEdf = melt(data.frame(UE = UE,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  EEUdf = melt(data.frame(EEU = EEU,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  U2Edf = melt(data.frame(U2E = U2E,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  U2E2df = melt(data.frame(U2E2 = U2E2,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  
  g1 = ggplot(EEdf,aes(x=lag,y=value,group=variable)) + geom_line(aes(linetype=variable)) + theme(legend.position="none",axis.title.x=element_blank()) +
    ylab(TeX('$\\phi_{\\xi\\xi}(\\tau)$')) + ylim(-1,1) + xlim(-maxlag,maxlag)
  g2 = ggplot(UEdf,aes(x=lag,y=value,group=variable)) + geom_line(aes(linetype=variable)) + theme(legend.position="none",axis.title.x=element_blank()) +
    ylab(TeX('$\\phi_{u\\xi}(\\tau)$')) + ylim(-1,1) + xlim(-maxlag,maxlag)
  g3 = ggplot(EEUdf,aes(x=lag,y=value,group=variable)) + geom_line(aes(linetype=variable)) + theme(legend.position="none",axis.title.x=element_blank()) +
    ylab(TeX('$\\phi_{\\xi(\\xi u)}(\\tau)$')) + ylim(-1,1) + xlim(0,maxlag)
  g4 = ggplot(U2Edf,aes(x=lag,y=value,group=variable)) + geom_line(aes(linetype=variable)) + theme(legend.position="none",axis.title.x=element_blank()) +
    ylab(TeX('$\\phi_{(u^2)\\prime \\xi}(\\tau)$')) + ylim(-1,1) + xlim(-maxlag,maxlag)
  g5 = ggplot(U2E2df,aes(x=lag,y=value,group=variable)) + geom_line(aes(linetype=variable)) + theme(legend.position="none") +
    ylab(TeX('$\\phi_{(u^2)\\prime \\xi^2}(\\tau)$')) + xlab(TeX('$\\tau')) + ylim(-1,1) + xlim(-maxlag,maxlag) 
  
  multiplot(g1,g2,g3,g4,g5,cols=1)
  
}