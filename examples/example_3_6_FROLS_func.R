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

out = frols(P,Y,rho)

###########################################################################
# FROLS END
###########################################################################

# print important information ---------------------------------------------

print("Selected terms")
print(colnames(out$Psel))
print("Estimated parameters")
print(out$th)
print(out$ERR*100)

