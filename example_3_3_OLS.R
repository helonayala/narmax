# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

# example 3.3- billings 2013

rm(list=ls()) # remove all vars and functions
while(!is.null(dev.list())) dev.off()     # clear all graphs
cat("\014")   # clean console

source("library_sysid.R")

# table 3.1
Mat = matrix(c(9,-5,5,-1.53,9.08,
               1,-1,8,-0.39,7.87,
               2,-5,6,-3.26,3.01,
               8,-2,0,0.36,5.98,
               0,0,9,0.13,9.05), # the data elements 
             nrow=5,              # number of rows 
             byrow = TRUE)        # fill matrix by rows

P = Mat[,1:4]
Y = Mat[,5]

# ordinary least squares solution
th_ls = ginv(P) %*% Y # see text in example 3.3, below eq. 3.22

# orthogonal least squares

# -- remove the predictors here
#P = P[,c(1,2,3)] # true theta

niter = ncol(P)

# out = CGS(P) # classical Gram-Schmidt
out = MGS(P) # modified Gram-Schmidt

W = out$Q
A = out$A

Alpha = t(W) %*% W

g = rep(0,niter)
for (i in 1:niter) {
  g[i] = (Y %*% W[,i]) / (W[,i] %*% W[,i])
}

g2 = solve(Alpha) %*% t(W) %*% Y

ERR = rep(0,niter)
for (i in 1:niter){
  ERR[i] = ( (Y %*% W[,i])^2 )/ ( (Y %*% Y) * (W[,i] %*% W[,i]) )
}

th_OLS = solve(A,g)

disp('OLS estimated parameters',th_OLS)
disp('ERR estimated parameters',ERR)
# CONCLUSION
# 1- the ERR values are not the same (confronting with table 3.2)
# 2- the values for th are the same for Table 3.2, varying the terms
# 3- it should be fine







