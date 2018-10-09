#' @title Plot correlation based tests
#'
#' @description See Billings book, chapter 5
#' @export
xcorrel = function(e,u) {

  maxlag = 25
  N = length(u)

  conf_factor = 1.96/sqrt(N)
  lag_vec = -maxlag:maxlag

  EE = crossco(e,e,maxlag)
  UE = crossco(u,e,maxlag)
  EEU = crossco(e[1:(N-1)],e[2:N]*u[2:N],maxlag)
  U2E = crossco(u^2 - mean(u^2),e,maxlag)
  U2E2 = crossco(u^2 - mean(u^2),e^2,maxlag)

  EEdf = reshape2::melt(data.frame(EE = EE,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  UEdf = reshape2::melt(data.frame(UE = UE,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  EEUdf = reshape2::melt(data.frame(EEU = EEU,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  U2Edf = reshape2::melt(data.frame(U2E = U2E,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")
  U2E2df = reshape2::melt(data.frame(U2E2 = U2E2,lag = lag_vec,conf_factor1 = conf_factor,conf_factor2 = -conf_factor),id.vars = "lag")

  g = list()

  g$g1 = ggplot2::ggplot(EEdf,ggplot2::aes(x=lag,y=value,group=variable)) + ggplot2::geom_line(ggplot2::aes(linetype=variable)) + ggplot2::theme(legend.position="none",axis.title.x=ggplot2:: element_blank()) +
    ggplot2::ylab(latex2exp::TeX('$\\phi_{\\xi\\xi}(\\tau)$')) + ggplot2:: ylim(-1,1) + ggplot2::xlim(-maxlag,maxlag)
  g$g2 = ggplot2::ggplot(UEdf,ggplot2::aes(x=lag,y=value,group=variable)) + ggplot2::geom_line(ggplot2::aes(linetype=variable)) + ggplot2::theme(legend.position="none",axis.title.x=ggplot2:: element_blank()) +
    ggplot2::ylab(latex2exp::TeX('$\\phi_{u\\xi}(\\tau)$')) + ggplot2:: ylim(-1,1) + ggplot2::xlim(-maxlag,maxlag)
  g$g3 = ggplot2::ggplot(EEUdf,ggplot2::aes(x=lag,y=value,group=variable)) + ggplot2::geom_line(ggplot2::aes(linetype=variable)) + ggplot2::theme(legend.position="none",axis.title.x=ggplot2:: element_blank()) +
    ggplot2::ylab(latex2exp::TeX('$\\phi_{\\xi(\\xi u)}(\\tau)$')) + ggplot2:: ylim(-1,1) + ggplot2::xlim(0,maxlag)
  g$g4 = ggplot2::ggplot(U2Edf,ggplot2::aes(x=lag,y=value,group=variable)) + ggplot2::geom_line(ggplot2::aes(linetype=variable)) + ggplot2::theme(legend.position="none",axis.title.x=ggplot2:: element_blank()) +
    ggplot2::ylab(latex2exp::TeX('$\\phi_{(u^2)\\prime \\xi}(\\tau)$')) + ggplot2:: ylim(-1,1) + ggplot2::xlim(-maxlag,maxlag)
  g$g5 = ggplot2::ggplot(U2E2df,ggplot2::aes(x=lag,y=value,group=variable)) + ggplot2::geom_line(ggplot2::aes(linetype=variable)) + ggplot2::theme(legend.position="none") +
    ggplot2::ylab(latex2exp::TeX('$\\phi_{(u^2)\\prime \\xi^2}(\\tau)$')) + ggplot2::xlab(latex2exp::TeX('$\\tau')) + ggplot2:: ylim(-1,1) + ggplot2::xlim(-maxlag,maxlag)

  return(g)
}
