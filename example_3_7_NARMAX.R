# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

library(tidyverse)
library(MASS)
source("library_sysid.R")
set.seed(0) # reproducibility

# code begins

# generate simulation data ------------------------------------------------
N = 400 # how many data?
u = rnorm(N,mean = 0,sd=1)
e = rnorm(N,mean = 0,sd=0.04^2)
y = rep(0,length(u))
for (k in 3:N) {
  y[k] = 0.5*y[k-1] + 
    u[k-2] +
    0.1*(u[k-2]^2) +
    0.5*e[k-1] +
    0.1*u[k-1]*e[k-2] +
    e[k]
}
e = rep(0,N)

# Step 1: define model initial conditions
# model parameters --------------------------------------------------------
rho_p = 1e-2
rho_n = 1.9e-6
#rho_n = 4.55e-7
nu = 2
ny = 2
ne = 2
l = 2
n = nu + ny + ne
p = max(nu,ny,ne)+1

Y = targetVec(y,p)

selectTerms = NULL # terms selected for parsimonious model from the NARMAX full-model 

# Step 2: identify process sub-model (NARX) FROLS
out = regMatNARMAX(u,y,e,nu,ny,ne,p,l,selectTerms)
Pp  = out$Pp

outNARX = frols(Pp,Y,rho_p)
e = c(rep(0,p-1), Y - outNARX$W %*% outNARX$g)

print("---------------------------------")
print("-- Step 2: NARX TERM SELECTION --")
print("---------------------------------")
print((outNARX$ERR)*100)
print("SUM ERR: ")
print(sum(outNARX$ERR)*100)
print("selected terms")
print(colnames(outNARX$Psel))

# Step 3: identify noise-related sub-model (MA) FROLS
out = regMatNARMAX(u,y,e,nu,ny,ne,p,l,selectTerms)
Pnp = cbind(outNARX$Psel,out$Pnp)
outNARMAX = frols(Pnp,Y,rho_n)

selectTerms = colnames(outNARMAX$Psel) # final model terms selection

print("-----------------------------------")
print("-- Step 3: NARMAX TERM SELECTION --")
print("-----------------------------------")
print((outNARMAX$ERR)*100)
print("SUM ERR: ")
print(sum(outNARMAX$ERR)*100)
print("selected terms")
print(selectTerms)
print("estimated parameters")
print(outNARMAX$th)

nth = length(outNARMAX$th)
e = c(rep(0,p-1), Y - outNARMAX$W %*% outNARMAX$g)

# Step 4: obtain final model parameters (NARMAX) ELS
iterELS = 10

Th_NARMAX_hat = matrix(0,nrow = nth,ncol = iterELS)
dlt = rep(0,iterELS)
for(s in 1:iterELS){
  
  out = regMatNARMAX(u,y,e,nu,ny,ne,p,l,selectTerms)
  P = out$P
  
  th_NARMAX_hat = ginv(P) %*% Y

  # calculate error and pad zeros for the initial conditions
  e_1 = e
  e = c(rep(0,p-1),Y - (P %*% th_NARMAX_hat)[,])
  
  # --- stop conditions
  dlt[s] = sqrt(sum((e - e_1)^2))
  
  # save estimated vectors
  Th_NARMAX_hat[,s] = th_NARMAX_hat
}
plot(dlt)
# final print

print("-----------------------------------")
print("-----------------------------------")
print("-----------------------------------")
print("End of NARMAX system identification procedure")
print("selected terms:")
print(colnames(P))
print("related parameters:")
print(Th_NARMAX_hat[,iterELS])







