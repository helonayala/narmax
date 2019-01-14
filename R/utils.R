
#' @title Subset matrix
#' @description Subset a matrix like mat[rows, cols]. Return a matrix even when the
#' subset generates a lower dimensional structure, i.e., never loose rownames or colnames.
#' @param mat Target matrix
#' @param rows Rows subsetter
#' @param cols Cols subsetter
#' @export
subsetMatrix = function (mat, rows, cols) {
  names = dimnames(mat)
  rows = if (is.null(rows)) 1:nrow(mat) else rows
  cols = if (is.null(cols)) 1:ncol(mat) else cols
  rownames = if (typeof(rows) == 'character') rows else names[[1]][rows]
  colnames = if (typeof(cols) == 'character') cols else names[[2]][cols]
  submat = matrix(
    mat[rows, cols],
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(rownames, colnames)
  )
  return(submat)
}

#' @title Find term indexes
#' @description Find term indexes in regression matrix
#' @param P Regression matrix
#' @param terms Vector of terms (string format)
findTermIndexes = function (P, terms) {
  names = colnames(P)
  indexes = which(names == terms[1])
  if (length(terms) > 1) {
    for (i in 2:length(terms)) {
      indexes = c(indexes, which(names == terms[i]))
    }
  }
  return(indexes)
}

#' @title Calculate normalized cross-correlation of 2 signals
#' Constructed on the basis of the code provided by Norgaard, Ravn, Poulsen and Hansen: https://www.mathworks.com/matlabcentral/fileexchange/87-nnsysid
#' @description See Billings book, chapter 5
#' @export
crossco = function(v,w,maxlag){
  v = v-mean(v)
  w = w-mean(w)

  nvw = length(v)

  normcoef = sqrt(sum(v*v)*sum(w*w))
  coefs = rep(0,maxlag+1)
  for (k in 0:maxlag){
    coefs[k+1] = sum(v[1:(nvw-k)]*w[(k+1):nvw])/normcoef
  }

  coefs2 = rep(0,maxlag)
  for(k in 1:maxlag){
    coefs2[k] = sum(w[1:(nvw-k)]*v[(k+1):nvw])/normcoef
  }
  coefs = c(rev(coefs2),coefs)

  return(coefs)
}


#' @title Calculate R2
#'
#' @description Calculate R2 (multiple correlation coefficient) for a predicted value
#' @param real measured data
#' @param est predicted data
#' @return R2
#' @export
calcR2 = function(real,est){
  SSE = sum((real - est)^2)
  avg_real = mean(real)
  sum2 = sum((real - avg_real)^2)
  R2 = 1 - SSE / sum2
  return(R2)
}

#' @title Prepare plots for predict() output
#'
#' @description Prepare function calls for output
#' @export
cookPlots <- function(df,R2,type) {

  df2 = tidyr::gather(df,variable, measurement, -time)

  dfy = dplyr::filter(df2,variable %in% c("y","yh"))

  dfy$label = ifelse(dfy$variable == "y", "Real", "Predicted")

  p1 = ggplot2::ggplot(data = dfy) +
    ggplot2::geom_line(ggplot2::aes(x=time, linetype=label,y = measurement)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.3),
      axis.title.x = ggplot2::element_text(vjust = - 0.5),
      plot.title = ggplot2::element_text(vjust = 1.5),
      legend.title = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      #legend.key.size = unit(0.5, "in"),
      legend.position = "bottom") +
    ggplot2::ggtitle(paste0('Predictions in ',type,". R2 = ",sprintf("%0.4f",R2)))+
    ggplot2::xlab("Time") + ggplot2::ylab("Output")

  p2 = ggplot2::ggplot(data = df, ggplot2::aes(x = time)) +
    ggplot2::geom_line(ggplot2::aes(y = e)) +
    ggplot2::scale_color_manual(values = c("#000000")) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste0('Error in ',type,". R2 = ",sprintf("%0.4f",R2)))+
    ggplot2::xlab("Time") +
    ggplot2::ylab("Error")

  out = list(p1 = p1,
             p2 = p2)

  return(out)
}

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
         { # 3 undefined
           df = NULL
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

#' @title Plot correlation based tests
#'
#' @description See Billings book, chapter 5
#' @export
xcorrel.ts = function(e) {

  maxlag = 25
  N = length(e)

  conf_factor = 1.96/sqrt(N)
  lag_vec = -maxlag:maxlag

  Ep  = e - mean(e)
  E2p = e^2 - mean(e^2)

  EpEp   = crossco(Ep,Ep,  maxlag)
  EpE2p  = crossco(Ep,E2p, maxlag)
  E2pE2p = crossco(E2p,E2p,maxlag)

  df = tidyr::gather(data.frame(lag = -maxlag:maxlag,EpEp,EpE2p,E2pE2p),
                     variable, measurement, -lag,factor_key = TRUE)

  levels(df$variable) =  c(latex2exp::TeX('$\\phi_{\\xi\\prime\\xi\\prime}(\\tau)$'),
                           latex2exp::TeX('$\\phi_{\\xi\\prime(\\xi^2)\\prime}(\\tau)$'),
                           latex2exp::TeX('$\\phi_{(\\xi^2)\\prime(\\xi^2)\\prime}(\\tau)$'))

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





