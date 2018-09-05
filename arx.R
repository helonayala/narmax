# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

# load libraries
library(ggplot2) # fancy plots
library(signal)  # filter for input signal
library(MASS)    # use ginv
library(tidyr)   # melt data frames
library(dplyr)   # create tibbles

# load functions
source("library_sysid.R")

# allows reproducibility
set.seed(42)
#model parameters
na = 4
nb = 3
p  = 1 + max(na,nb)
sd_noise = 1e-1

# generate input-output data ----------------------------------------------
N = 1000 # number of samples
cutoff = 0.1 # normalized frequency cutoff

# create input signal
ue = multisine(N,cutoff)
uv = randnoise(N,cutoff)

# M_spec(ue,'ue')
# M_spec(uv,'uv')

ye = array(0,N)
yv = array(0,N)

for (k in p:N){
  ye[k] = -0.3*ye[k-1] -0.5*ye[k-2] -0.1*ye[k-3] - 0.4*ye[k-4] + 0.2*ue[k-1] + 0.32*ue[k-2] + 0.1*ue[k-3]
  yv[k] = -0.3*yv[k-1] -0.5*yv[k-2] -0.1*yv[k-3] - 0.4*yv[k-4] + 0.2*uv[k-1] + 0.32*uv[k-2] + 0.1*uv[k-3]
}
yeor = ye
ye = rnorm(N,mean=ye,sd=sd_noise)
yvor = yv
yv = rnorm(N,mean=yv,sd=sd_noise)

Phie = regMatrix(ye,ue,na,nb)
Ye   = targetVec(ye,na,nb)
Phiv = regMatrix(yv,uv,na,nb)
Yv   = targetVec(yv,na,nb)

# estimate parameters -----------------------------------------------------
th_hat = ginv(Phie) %*% Ye

# calculate predictions ---------------------------------------------------
ye_osa = (Phie %*% th_hat)[,]
ye_fr = calcFR_ARX(ye,ue,na,nb)
yv_osa = (Phiv %*% th_hat)[,]
yv_fr = calcFR_ARX(yv,uv,na,nb)

# create tibbles for plotting ---------------------------------------------

df_dataset = tibble(time = 1:N,
                    ye = ye,
                    yeor = yeor,
                    ue = ue,
                    yv = yv,
                    yvor = yvor,
                    uv = uv) %>%
  gather(variable, measurement, -time)

df_pred = tibble(time = p:N,
                 ye_osa = ye_osa,
                 ye_fr  = ye_fr,
                 yv_osa = yv_osa,
                 yv_fr  = yv_fr) %>%
  gather(variable, measurement, -time)

df_error = tibble(time = p:N,
                  ee_osa = Ye - ye_osa,
                  ee_fr  = Ye - ye_fr,
                  ev_osa = Yv - yv_osa,
                  ev_fr  = Yv - yv_fr) %>%
  gather(variable, measurement, -time)

df_all = bind_rows(df_dataset,df_pred)


# plots -------------------------------------------------------------------
p = list()

p = c(p,list(
  ggplot(data=filter(df_dataset,variable %in% c("ue","uv"))) +
    geom_line(aes(x = time,y =measurement,color=variable)) +
    ggtitle("Input"))
)

p = c(p,list(
  ggplot(data=filter(df_dataset,variable %in% c("ye","yeor"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("Output (estimation)"))
)

p = c(p,list(
  ggplot(data=filter(df_dataset,variable %in% c("yv","yvor"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("Output (valitation)"))
)

p = c(p,list(
  ggplot(data=filter(df_all,variable %in% c("ye","ye_osa","ye_fr"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("predictions (estimation"))
)

p = c(p,list(
  ggplot(data=filter(df_all,variable %in% c("yv","yv_osa","yv_fr"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("predictions (validation"))
)

p = c(p,list(
  ggplot(data=filter(df_error,variable %in% c("ee_osa","ee_fr"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("prediction residuals (estimation)"))
)

p = c(p,list(
  ggplot(data=filter(df_error,variable %in% c("ev_osa","ev_fr"))) + 
    geom_line(aes(x = time,y =measurement,color=variable)) + 
    ggtitle("prediction residuals (validation)"))
)

print(p)






