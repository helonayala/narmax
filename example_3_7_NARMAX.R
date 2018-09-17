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

# Step 1: define model initial conditions
# model parameters --------------------------------------------------------
rho_p = 0.05
rho_n = 0.01
nu = 2
ny = 2
ne = 2
l = 2
n = nu + ny + ne
p = max(nu,ny,ne)+1

plot(y)
lines(y)

# Step 2: identify process sub-model (NARX)



# Step 3: identify noise-related sub-model (MA)



# Step 4: obtain final model parameters (ELS)

