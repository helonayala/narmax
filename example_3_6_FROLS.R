# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

library(tidyverse)
source("library_sysid.R")
set.seed(0) # reproducibility

# code begins

# generate simulation data ------------------------------------------------
N = 200
u = runif(N,min = -1,max = 1)
e = rnorm(N,mean = 0,sd=0.1)
y = rep(0,length(u))
for (k in 3:N) {
  y[k] = -0.605*y[k-1] - 0.163*y[k-2]^2 + 0.588*u[k-1] - 0.24*u[k-2] + e[k]
}

# model parameters --------------------------------------------------------
rho = 0.04
nu = 2
ny = 2
ne = 0
l = 3
n = nu + ny + ne
p = max(nu,ny,ne)+1

# create regression and target matrices -----------------------------------
P = regMatNARX(u,y,nu,ny,l)
Y = targetVec(y,p)

M = ncol(P)
NP = nrow(P)

###########################################################################
# FROLS BEGIN
###########################################################################

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
  
  Qs = cbind(Qs,Qm[,ls])
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
###########################################################################
# FROLS END
###########################################################################


# print important information ---------------------------------------------

print("Selected terms")
print(colnames(Psel))
print("Estimated parameters")
print(th_FROLS[,])
print(ERRvec[,]*100)

