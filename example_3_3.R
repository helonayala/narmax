# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

# example 3.3

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

# table 3.1
Mat = matrix(c(9,-5,5,-1.53,9.08,
               1,-1,8,-0.39,7.87,
               2,-5,6,-3.26,3.01,
               8,-2,0,0.36,5.98,
               0,0,9,0.13,9.05), # the data elements 
             nrow=5,              # number of rows 
             byrow = TRUE)        # fill matrix by rows

Phi = Mat[,1:4]
Y  = Mat[,5]

# ordinary least squares solution
th_ls = ginv(Phi) %*% Y # see text in example 3.3, below eq. 3.22

# orthogonal least squares

