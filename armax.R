# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls())
cat("\014")  

# load libraries
library(ggplot2) # fancy plots
library(signal)  # filter for input signal
library(MASS)

# load functions
source("library_sysid.R")

# allows reproducibility
set.seed(42)
#model parameters
N = 1000 # number of samples
na = 4
nb = 3
p  = 1 + max(na,nb)
sd_noise = 1e-2
# create input signal
bf = butter(6, 1/10, type="low")
u  = filter(bf, rnorm(N,mean=0,sd=1))

plot(u)

y = array(0,N)

for (k in p:N){
  y[k] = -0.3*y[k-1] -0.5*y[k-2] -0.1*y[k-3] - 0.4*y[k-4] + 0.2*u[k-1] + 0.32*u[k-2] + 0.1*u[k-3]
}
yor = y
y = rnorm(N,mean=y,sd=sd_noise)

Phi = regMatrix(y,u,na,nb)
Y   = targetVec(y,na,nb)

th_hat = ginv(Phi) %*% Y










