# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls())
cat("\014")  

library(ggplot2)
source("library_sysid.R")

y = 1:10
u = 1:10

Phi = regMatrix(y,u,2,3)

