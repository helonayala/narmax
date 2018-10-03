#' @title Clear Workspace
#'
#' @description Clear the workspace (Actually just a rdoxygen example)
#' @export
clearWorkspace = function () {
  rm(list=ls())
  cat('\014')
  while (!is.null(dev.list())) dev.off()
}

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

#' @title Decibel
#'
#' @description See https://en.wikipedia.org/wiki/Decibel
#' @param X input
#' @return db, X in db
#' @export
db = function(X) {

  db = 20*log10(abs(X))

  return(db)
}

#' @title Generate filtered random noise signal excitation signal
#'
#' @description Filters N-normally distributed samples with a defined cut-off frequency
#'
#' Ref.: Schoukens, Johan, Rik Pintelon, and Yves Rolain. Mastering system identification in 100 exercises. John Wiley & Sons, 2012.
#' @param N number of samples
#' @param cutoff cutoff normalized frequency
#' @return u, excitation signal
#' @export
randnoise = function(N,cutoff){

  bf = signal::butter(6, cutoff, type="low")
  u  = signal::filter(bf, rnorm(N,mean=0,sd=1))

  return (u)
}

#' @title Generate multisine excitation signal
#'
#' @description Ref.: Schoukens, Johan, Rik Pintelon, and Yves Rolain. Mastering system identification in 100 exercises. John Wiley & Sons, 2012.
#' @param N number of samples
#' @param cutoff cutoff frequency (normalized)
#' @return u, excitation signal
#' @export
multisine = function(N,cuoff){

  # number of sines
  Ns = round(N*cutoff)
  # normalized frequencies
  f = (0:(N-1))/N

  # create ifft argument
  U = matrix(0,N,1)
  U[2:(Ns+1)] = exp(1i*2*pi*runif(Ns))

  # perform ifft
  u = 2*Re(signal::ifft(U))
  u = u/sd(u)

  return(u[,])
}

#' @title ploting signal spectrum
#'
#' @description performs U = fft(u)/sqrt(N), and plots 20*log10(abs(U)) for f in [0,0.5] (normalized)
#' @param u input vector to be analyzed
#' @examples
#' t=seq(from = 0, to = 400*2*pi, by = 0.01)
#' t = t[1:2^14]
#' u=sin(2*pi*t)
#' M_spec(u)
#' @export
M_spec = function(u){

  N = length(u)

  U = stats::fft(u)

  U = U[1:round(N/2)]/sqrt(N)

  f = ( (1:round(N/2)) - 1)/N

  plot(f,db(U),xlab =  'Frequency (normalized)', ylab =  'Amplitude (dB)')
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

#' @title Multiple plot function
#'
#' @description Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2::ggplot2)/
#' ggplot2::ggplot objects can be passed in ..., or to plotlist (as a list of ggplot2::ggplot objects)
#' - cols:   Number of columns in layout
#' - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


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

#' @title Calculate normalized cross-correlation of 2 signals
#'
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
