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
  # return(coefs)
}