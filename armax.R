# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls())
cat("\014")  

# load libraries
library(ggplot2) # fancy plots
library(signal)  # filter for input signal
library(MASS)
library(tidyr)

# load functions
source("library_sysid.R")

# allows reproducibility
set.seed(42)
#model parameters
N = 1000 # number of samples
na = 4
nb = 3
p  = 1 + max(na,nb)
sd_noise = 1e-1
# create input signal
bf = butter(6, 1/10, type="low")
u  = signal::filter(bf, rnorm(N,mean=0,sd=1))

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

# calculate predictions
y_osa = (Phi %*% th_hat)[,]
y_fr = y[1:(p-1)]
u_fr = u[1:(p-1)]
for (k in p:N){
  phi_k = regMatrix(c(y_fr[(k-p+1):(k-1)],0),c(u_fr[(k-p+1):(k-1)],0),na,nb)
  y_fr[k] = phi_k %*% th_hat
  u_fr[k] = u[k]
}
y_fr = y_fr[p:N]

df_dataset = tibble(time = 1:N,
                        y = y,
                        yor = yor,
                        u = u) %>%
  gather(variable, measurement, -time)

df_pred = tibble(time = p:N,
                 y_osa = y_osa,
                 y_fr = y_fr) %>%
  gather(variable, measurement, -time)

df_error = tibble(time = p:N,
                 e_osa = Y - y_osa,
                 e_fr = Y - y_fr) %>%
  gather(variable, measurement, -time)


df_all = bind_rows(df_dataset,df_pred)

# ggplot(data=filter(df_dataset,variable %in% c("y","yor"))) + 
#   geom_line(aes(x = time,y =measurement,color=variable)) + 
#   ggtitle("Output")

ggplot(data=filter(df_all,variable %in% c("yor","y_osa","y_fr"))) + 
  geom_line(aes(x = time,y =measurement,color=variable)) + 
  ggtitle("predictions")

ggplot(data=filter(df_error,variable %in% c("e_osa","e_fr"))) + 
  geom_line(aes(x = time,y =measurement,color=variable)) + 
  ggtitle("prediction residuals")

# ggplot(data=filter(df_dataset,variable %in% c("u"))) + 
#   geom_line(aes(x = time,y =measurement,color=variable)) + 
#   ggtitle("Input")









