#' @title Plot correlation based tests
#'
#' @description See Billings book, chapter 5
#' @export
xcorrel = function(e,u,nl) {

  maxlag = 25
  N = length(u)

  conf_factor = 1.96/sqrt(N)
  lag_vec = -maxlag:maxlag

  EE = crossco(e,e,maxlag)
  UE = crossco(u,e,maxlag)
  EEU = crossco(e[1:(N-1)],e[2:N]*u[2:N],maxlag)
  U2E = crossco(u^2 - mean(u^2),e,maxlag)
  U2E2 = crossco(u^2 - mean(u^2),e^2,maxlag)

  df1 = tidyr::gather(data.frame(lag = -maxlag:maxlag,EE,UE),
                      variable, measurement, -lag,factor_key = TRUE)
  df2 = tidyr::gather(data.frame(lag = 0:maxlag, EEU = EEU[(maxlag+1):(2*maxlag+1)]),
                      variable, measurement, -lag,factor_key = TRUE)
  df3 = tidyr::gather(data.frame(lag = -maxlag:maxlag,U2E,U2E2),
                      variable, measurement, -lag,factor_key = TRUE)

  if (class(mdl) %in% c("armax","arx")){
    nl = 0
  } else  if (class(mdl) %in% "narmax"){
    nl = 1
  } else {
    nl = 2
  }
  switch(nl,
    { # 1 linear with X
      df = rbind(df1)
      levels(df$variable) =  c(latex2exp::TeX('$\\phi_{\\xi\\xi}(\\tau)$'),
                               latex2exp::TeX('$\\phi_{u\\xi}(\\tau)$'))
    },
    { # 2 nonlinear with x
      df = rbind(df1,df2,df3)
      levels(df$variable) =  c(latex2exp::TeX('$\\phi_{\\xi\\xi}(\\tau)$'),
                               latex2exp::TeX('$\\phi_{u\\xi}(\\tau)$'),
                               latex2exp::TeX('$\\phi_{\\xi(\\xi u)}(\\tau)$'),
                               latex2exp::TeX('$\\phi_{(u^2)\\prime \\xi}(\\tau)$'),
                               latex2exp::TeX('$\\phi_{(u^2)\\prime \\xi^2}(\\tau)$'))
    },
    { # 3 nonlinear ts
      df = data.frame() # WIP
    })

  g = ggplot2::ggplot(df) +
    ggplot2::geom_line(ggplot2::aes(x = lag, y = measurement)) + #,color = variable)) +
    ggplot2::geom_hline(yintercept=conf_factor, linetype="dashed") + ggplot2::geom_hline(yintercept=-conf_factor, linetype="dashed") +
    ggplot2::ylim(-1,1)  +
    ggplot2::ylab("") +
    ggplot2::facet_grid(~variable, scales = "free",labeller = ggplot2::label_parsed) +
    #ggplot2::theme(strip.text.x = ggplot2::element_text(size = 20))+
    ggplot2::theme_bw(base_size = 14)

  return(g)
}
