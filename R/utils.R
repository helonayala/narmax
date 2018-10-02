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
  if (is.null(rows)) rows = 1:nrow(mat)
  if (is.null(cols)) cols = 1:ncol(mat)
  names = dimnames(mat)
  submat = matrix(
    mat[rows, cols],
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(names[[1]][rows], names[[2]][cols])
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
  for (i in 2:length(terms)) {
    indexes = c(indexes, which(names == terms[i]))
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
#' @examples
#' t=seq(from = 0, to = 400*2*pi, by = 0.01)
#' t = t[1:2^14]
#' u=sin(2*pi*t)
#' M_spec(u)
#' @export
calcR2 = function(real,est){

  SSE = sum((real - est)^2)

  avg_real = mean(real)

  sum2 = sum((real - avg_real)^2)

  R2 = 1 - SSE / sum2

  return(R2)
}

