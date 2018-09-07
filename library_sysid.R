regMatrix_ARX = function(y,u,na,nb,p){
  # creates the regression matrix
  
  N = length(y)
  
  if (N != length(u)) {
    stop("input-output vectors should have the same size")
  }

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

regMatrix_MA = function(y,u,e,na,nb,nc,p){
  # creates the regression matrix
  
  N = length(y)
  
  if (N != length(u) || N != length(e))  {
    stop("input vectors should have the same size")
  }
  
  Phi = matrix(0, nrow = N-p+1, ncol = na+nb+nc)
  
  colPhi = NULL
  for(i in 1:na){
    Phi[,i] = -y[(p-i):(N-i)]
    colPhi = c(colPhi,paste0("-y(k-",i,")"))
  }
  for(i in 1:nb){
    Phi[,na+i] = u[(p-i):(N-i)]
    colPhi = c(colPhi,paste0("u(k-",i,")"))
  }
  for(i in 1:nc){
    Phi[,na+nb+i] = e[(p-i):(N-i)]
    colPhi = c(colPhi,paste0("e(k-",i,")"))
  }
  
  rowPhi = paste0(rep("k=",N-p+1),p:N)
  
  colnames(Phi) = colPhi
  rownames(Phi) = rowPhi
  
  return(Phi)
}

targetVec = function(y,p){
  # returns the target vector
  
  N = length(y)
    
  Y = y[p:N]
  return(Y)
}

calcFR_ARX = function(y,u,na,nb,p,th_hat){
  
  y_fr = y[1:(p-1)]
  u_fr = u[1:(p-1)]
  
  for (k in p:N){
    phi_k = regMatrix_ARX(c(y_fr[(k-p+1):(k-1)],0),c(u_fr[(k-p+1):(k-1)],0),na,nb,p)
    y_fr[k] = phi_k %*% th_hat
    u_fr[k] = u[k]
  }
  y_fr = y_fr[p:N]
  
  return(y_fr)
}

calcOSA_ARMAX = function(y,u,na,nb,nc,p,th_hat){
  
  y_fr = y[1:(p-1)]
  u_fr = u[1:(p-1)]
  e_fr = rep(0,p-1)
  
  for (k in p:N){
    auxy = c(y_fr[(k-p+1):(k-1)],0)
    auxu = c(u_fr[(k-p+1):(k-1)],0)
    auxe = c(e_fr[(k-p+1):(k-1)],0)
    phi_k = regMatrix_MA(auxy,auxu,auxe,na,nb,nc,p)
    y_fr[k] = phi_k %*% th_hat
    u_fr[k] = u[k]
    e_fr[k] = y[k] - y_fr[k]
  }
  y_fr = y_fr[p:N]
  
  return(y_fr)
}

db = function(X){
  X = abs(X)^2
  # from matlab function:
  # We want to guarantee that the result is an integer
  # if X is a negative power of 10.  To do so, we force
  # some rounding of precision by adding 300-300.
  
  Y = (10*log10(X)+300)-300;
  
}

randnoise = function(N,cutoff){
  # generate low-pass filtered random noise 
  
  # N number of samples
  # cutoff cutoff normalized frequency
  
  # adapted from Pintelon,Schoukens book
  
  bf = butter(6, cutoff, type="low")
  u  = signal::filter(bf, rnorm(N,mean=0,sd=1))
  
  return (u)
}

multisine = function(Ndata1,cuoff){
  # generate low-pass filtered random noise 
  
  # N number of samples
  # cutoff cutoff normalized frequency
  
  # adapted from Pintelon,Schoukens book
  
  fsample = 1 # use normalized frequency
  Ts = 1/fsample
  Nsines = round(Ndata1*cutoff)
  f=(0:(Ndata1-1))*fsample/Ndata1
  
  U = matrix(0,Ndata1,1)
  U[2:(Nsines+1)] = exp(1i*2*pi*runif(Nsines));
  u = 2*Re(ifft(U))
  u = u/sd(u)

  return(u[,])
}

M_spec = function(u,title = 'u'){
  # returns matrix for ploting signal spectra
  
  # adapted from Pintelon,Schoukens book
  
  Ndata1 = length(u)
  
  Um=fft(u)/sqrt(Ndata1)
  
  f = (0:(Ndata1-1))/Ndata1
  LinesPlot=(1:floor(Ndata1/2))
  
  plot(f[LinesPlot],db(Um[LinesPlot]),xlab =  'Frequency (normalized)', ylab =  'Amplitude (dB)', main = title)
}

calcR2 = function(real,est){
  
  SSE = sum((real - est)^2)
  
  avg_real = mean(real)
  
  sum2 = sum((real - avg_real)^2)
  
  R2 = 1 - SSE / sum2
  
  return(R2)
}  
  

CGS = function(P) {
  # Classical Gram-Schimidt factorization
  # obtains P = Q * A
  # where P is N x Nth
  # so that
  # Q is a N x Nth matrix with orthogonal columns 
  # A is a Nth x Nth unit upper triangular matrix
  
  N   = nrow(P)
  Nth = ncol(P)
  
  # init matrix
  A = eye(Nth)
  Q = matrix(0,nrow =N, ncol = Nth)
  Q[,1] = P[,1]
  for (i in 2:Nth) {
    Q[,i] = P[,i] 
    for (j in 1:(i-1)){
      disp(j,i)
      A[j,i] = (Q[,j] %*% P[,i]) / (Q[,j] %*% Q[,j])
      Q[,i] = Q[,i] - A[j,i] * Q[,j]
    }
  }
  # THE END
  
  # format output
  out$A = A
  out$Q = Q
  return(out)
}
