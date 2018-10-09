
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
