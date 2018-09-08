# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

source("library_sysid.R")

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
nu = 2
ny = 3
ne = 0
l = 3
n = nu + ny + ne
p = max(nu,ny,nu)+1
# candidate
ncan = 0
canl = NULL

Phi = regMatNARX(u,y,nu,ny,l)
Y = targetVec(y,p)

