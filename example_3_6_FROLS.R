# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

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
ny = 2
ne = 0
l = 3
n = nu + ny + ne
# candidate
ncan = 0
canl = NULL

canl3 = expand.grid(1:n,1:n,1:n)
canl2 = expand.grid(1:n,1:n)
canl1 = expand.grid(1:n)

canl3_a = t(apply(canl3,1,sort)) # order each row
canl3_b = unique(canl3_a) # keep unique values
canl2_a = t(apply(canl2,1,sort)) # order each row
canl2_b = unique(canl2_a) # keep unique values
canl1_a = t(apply(canl1,1,sort)) # order each row
canl1_b = unique(canl1_a) # keep unique values


canl1_b
canl2_b
canl3_b