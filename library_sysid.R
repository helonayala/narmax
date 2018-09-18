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
    
  Y = matrix(y[p:N],ncol = 1)
  
  rownames(Y) = paste0(rep("k=",N-p+1),p:N)
  colnames(Y) = "y(k)"
  
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
  # Aguirre 2015 book
  # obtains P = Q * A
  # where P is N x Nth
  # so that
  # Q is a N x Nth matrix with orthogonal columns 
  # A is a Nth x Nth unit upper triangular matrix
  
  N   = nrow(P)
  Nth = ncol(P)
  
  # init matrix
  A = diag(1,nrow = Nth)
  Q = matrix(0,nrow =N, ncol = Nth)
  Q[,1] = P[,1]
  for (i in 2:Nth) {
    Q[,i] = P[,i] 
    for (j in 1:(i-1)){
      # disp(j,i)
      A[j,i] = (Q[,j] %*% P[,i]) / (Q[,j] %*% Q[,j])
      Q[,i] = Q[,i] - A[j,i] * Q[,j]
    }
  }
  # THE END
  
  out = list()
  # format output
  out$A = A
  out$Q = Q
  return(out)
}


MGS = function(P) {
  # Modified Gram-Schimidt factorization
  # Aguirre 2015 book
  # obtains P = Q * A
  # where P is N x Nth
  # so that
  # Q is a N x Nth matrix with orthogonal columns 
  # A is a Nth x Nth unit upper triangular matrix
  
  N   = nrow(P)
  Nth = ncol(P)
  
  # init matrix
  A = diag(1,nrow = Nth)
  P_i_1 = P
  Q = matrix(0,nrow =N, ncol = Nth)
  P_i  = matrix(0,nrow = N,ncol = Nth)
  for (i in 1:(Nth-1)) {
    Q[,i] = P_i_1[,i]
    for (j in (i+1):Nth){
      #disp(j,i)
      A[i,j] = (Q[,i] %*% P_i_1[,j]) / (Q[,i] %*% Q[,i])
      P_i[,j] = P_i_1[,j] - A[i,j] * Q[,i]
    }
    P_i_1 = P_i
  }
  # THE END
  Q[,Nth] = P_i_1[,Nth]
  
  out = list()
  # format output
  out$A = A
  out$Q = Q
  return(out)
}

regMatNARX = function(u,y,nu,ny,l){
  n = nu+ny
  p = max(nu,ny)+1
  auxexp = list()
  candlist = list()
  # generate all terms product combinations possible
  for(i in 1:l){
    eval(parse(text=paste0("auxexp$x",i,"=1:n"))) # generate input args for expand.grid
    cand = do.call(expand.grid,auxexp) # call expand.grid for changing number of arguments
    
    cand = t(apply(cand,1,sort)) # order each row of the matrix
    cand = unique(cand) # keep unique rows
    candlist[[i]] = cand
  }
  
  P0 = regMatrix_ARX(y,u,ny,nu,p) # creates the initial regression matrix (with ARX terms)
  P0[,1:ny] = -P0[,1:ny] # fix -y(k-..) from arx modeling
  colnames(P0) = c(paste0("y(k-",1:ny,")"),paste0("u(k-",1:nu,")")) # fix column names
  NP = nrow(P0)
  P = NULL # here we store the final NARX regression matrix 
  P = cbind(rep(1,NP), P) # constant term
  colnames(P) = "constant"
  P = cbind(P, P0) # l = 1, trivial
  
  for (i in 2:l){
    ncand = nrow(candlist[[i]])
    for(j in 1:ncand){
      Pcand_a = P0[,candlist[[i]][j,]]
      Pcand_b = matrix(apply(Pcand_a,1,prod),ncol = 1)
      colnames(Pcand_b) = str_c(colnames(P0[,candlist[[i]][j,]]),collapse = "")
      P = cbind(P,Pcand_b)
    }
  }
  return(P)
}

frols = function(P,Y,rho){
  
  
  ###########################################################################
  # FROLS BEGIN
  ###########################################################################
  
  M = ncol(P)
  NP = nrow(P)
  
  # 1st step ----------------------------------------------------------------
  sig = Y[,] %*% Y[,]
  selectTerms = NULL
  ERRvec = NULL
  gvec = NULL
  
  Qs = P
  g = rep(0,M)
  ERR = rep(0,M)
  for (m in 1:M){
    g[m] = (Y[,] %*% Qs[,m]) / (Qs[,m] %*% Qs[,m])
    ERR[m] = ( g[m]^2 * (Qs[,m] %*% Qs[,m]) ) / sig
  }
  l1 = which(ERR==max(ERR))
  selectTerms = l1 # vector keeping all selected terms
  
  # init
  # A = diag(1,M)
  A=1
  Qs = matrix(P[,l1],ncol = 1)
  gvec = g[l1]
  ERRvec = ERR[l1]
  
  # s-th step --------------------------------------------------------------
  for (s in 2:M){
    gm  = rep(0,M)
    Qm  = matrix(0,NP,M)
    ERR = rep(0,M)
    A = cbind(rbind(A,0),0)
    
    for (m in (1:M)[-selectTerms]) {
      
      sumQm = rep(0,NP)
      for (r in 1:(s-1)){
        sumQm = sumQm + ((P[,m] %*% Qs[,r]) /  (Qs[,r] %*% Qs[,r])) %*% Qs[,r]
      }
      Qm[,m] = P[,m] - sumQm
      
      gm[m] = (Y[,] %*% Qm[,m]) / (Qm[,m] %*% Qm[,m])
      
      ERR[m] = ( gm[m]^2 * (Qm[,m] %*% Qm[,m]) ) / sig
      
    }
    
    ls = which(ERR==max(ERR))
    selectTerms = cbind(selectTerms,ls) # vector keeping all selected terms
    
    Qs = cbind(Qs,Qm[,ls]) # keep set of orthogonal bases
    gvec = rbind(gvec,gm[ls])
    for (r in 1:(s-1)){
      A[r,s] = (Qs[,r] %*% P[,ls]) / (Qs[,r] %*% Qs[,r])
    }
    A[s,s] = 1
    
    ERRvec = rbind(ERRvec,ERR[ls])
    
    ESR = 1-sum(ERRvec)
    
    if (ESR <= rho){
      M0 = s
      break
    }
  }
  
  th_FROLS = solve(A,gvec)
  Psel = P[,selectTerms]
  
  out = list()
  out$th = th_FROLS[,]
  out$Psel = Psel
  out$g = gvec
  out$W = Qs
  out$A = A
  out$ERR = ERRvec[,]
  
  return(out)
  
}


regMatNARMAX = function(u,y,e,nu,ny,ne,p,l,selectTerms){
  n = nu+ny+ne
  P = NULL # here we store the final NARMAX regression matrix 
  
  P0 = regMatrix_MA(y,u,e,ny,nu,ne,p) # creates the initial regression matrix (with ARX terms)
  P0[,1:ny] = -P0[,1:ny] # fix -y(k-..) from armax modeling
  colnames(P0) = c(paste0("y(k-",1:ny,")"),paste0("u(k-",1:nu,")"),paste0("e(k-",1:ne,")")) # fix column names
  NP = nrow(P0)
  P = cbind(rep(1,NP), P) # constant term
  colnames(P) = "constant"
  P = cbind(P, P0) # l = 1, trivial
  
  # generate all terms product combinations possible
  auxexp = list()
  candlist = list()
  for(i in 1:l){
    eval(parse(text=paste0("auxexp$x",i,"=1:n"))) # generate input args for expand.grid
    cand = do.call(expand.grid,auxexp) # call expand.grid for changing number of arguments
    
    cand = t(apply(cand,1,sort)) # order each row of the matrix
    cand = unique(cand) # keep unique rows
    candlist[[i]] = cand
  }
  
  for (i in 2:l){
    ncand = nrow(candlist[[i]])
    for(j in 1:ncand){
      Pcand_a = P0[,candlist[[i]][j,]]
      Pcand_b = matrix(apply(Pcand_a,1,prod),ncol = 1)
      colnames(Pcand_b) = str_c(colnames(P0[,candlist[[i]][j,]]),collapse = "")
      P = cbind(P,Pcand_b)
    }
  }
  
  if (!is.null(selectTerms)){ # remove terms as indicated by FROLS (ELS part)
    ind2 = which(colnames(P) %in% selectTerms)
    P = P[,ind2]
  }
  
  ind = grepl("e(",colnames(P), fixed=TRUE) # find cols which have e[k-i] terms
  
  out = list()
  out$P = P
  out$Pp = P[,!ind]
  out$Pnp = P[,ind]
  
  return(out)
}







