regMatrix = function(y,u,na,nb){
  
  N = length(y)
  
  if (N != length(u)) {
    stop("input-output vectors should have the same size")
  }
  
  p = 1 + max(na,nb)
  
  Phi = matrix(0, nrow = N-p+1, ncol = na+nb)
  
  colPhi = NULL
  for(i in 1:na){
    Phi[,i] = -y[(p-i):(N-i)]
    colPhi = c(colPhi,paste0("-y(k-",i,")"))
  }
  for(i in 1:nb){
    Phi[,na+i] = u[(p-i):(N-i)]
    colPhi = c(colPhi,paste0("u(k-",i,")"))
  }
  
  rowPhi = paste0(rep("k=",N-p+1),p:N)
  
  colnames(Phi) = colPhi
  rownames(Phi) = rowPhi
  
  return(Phi)
}