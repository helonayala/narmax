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
nc = 2
p  = 1 + max(na,nb,nc)
sd_noise = 1e-2

# generate input-output data ----------------------------------------------
N = 1000 # number of samples
cutoff = 0.1 # normalized frequency cutoff
th = c(0.3,0.5,0.1,0.4,0.2,0.32,0.1)
# create input signal
ue = multisine(N,cutoff)
uv = randnoise(N,cutoff)

# M_spec(ue,'ue')
# M_spec(uv,'uv')

ye = array(0,N)
yv = array(0,N)

for (k in 5:N){
  phie = c(-ye[k-1],-ye[k-2],-ye[k-3],-ye[k-4],ue[k-1],ue[k-2],ue[k-3])
  phiv = c(-yv[k-1],-yv[k-2],-yv[k-3],-yv[k-4],uv[k-1],uv[k-2],uv[k-3])
  ye[k] = phie %*% th
  yv[k] = phiv %*% th
}
yeor = ye
ye = rnorm(N,mean=ye,sd=sd_noise)
yvor = yv
yv = rnorm(N,mean=yv,sd=sd_noise)

Phie = regMatrix_ARX(ye,ue,na,nb,p)
Ye   = targetVec(ye,p)

Phiv = regMatrix_ARX(yv,uv,na,nb,p)
Yv   = targetVec(yv,p)

# estimate parameters -----------------------------------------------------
niter = 10
Th_ARMAX_hat = matrix(0,na+nb+nc,niter)
th_ARX_hat = ginv(Phie) %*% Ye
th_ARX_hat0 = th_ARX_hat
th_ARMAX_hat0 = c(th_ARX_hat0,rep(0,nc))
ee_s1 = c(rep(0,p-1),Phie %*% th_ARX_hat)
dlt1 = rep(0,niter)
dlt2 = rep(0,niter)

for (i in 1:niter){
  
  Phie_ext = regMatrix_MA(ye,ue,ee_s1,na,nb,nc,p)
  
  th_ARMAX_hat = ginv(Phie_ext) %*% Ye
  
  # --- stop conditions
  dlt1[i] = sqrt(sum((th_ARMAX_hat - th_ARMAX_hat0)^2))
  th_ARMAX_hat0 = th_ARMAX_hat
  
  # calculate error and pad zeros for the initial conditions
  ee_s = ee_s1
  ee_s1 = c(rep(0,p-1),Ye - (Phie_ext %*% th_ARMAX_hat)[,])
  dlt2[i] = sqrt(sum((ee_s1 - ee_s)^2))
  
  # save estimated vectors
  Th_ARMAX_hat[,i] = th_ARMAX_hat
  th_ARX_hat = th_ARMAX_hat[1:(na+nb)]
}

plot(dlt1,main="armax theta convergence (log scale)",log="y")
plot(dlt2,main="armax error convergence (log scale)",log="y")

# calculate predictions ---------------------------------------------------
ye_osa = calcOSA_ARMAX(ye,ue,na,nb,nc,p,th_ARMAX_hat)
ye_fr  = calcFR_ARX(ye,ue,na,nb,p,th_ARX_hat)
yv_osa = calcOSA_ARMAX(yv,uv,na,nb,nc,p,th_ARMAX_hat)
yv_fr  = calcFR_ARX(yv,uv,na,nb,p,th_ARX_hat)

# prediction performance
R2e_osa = calcR2(Ye,ye_osa)
R2e_fr  = calcR2(Ye,ye_fr)
R2v_osa = calcR2(Yv,yv_osa)
R2v_fr  = calcR2(Yv,yv_fr)

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

# p = c(p,list(
#   ggplot(data=filter(df_dataset,variable %in% c("ue","uv"))) +
#     geom_line(aes(x = time,y =measurement,color=variable)) +
#     ggtitle("Input"))
#   )

# p = c(p,list(
#   ggplot(data=filter(df_dataset,variable %in% c("ye","yeor"))) + 
#   geom_line(aes(x = time,y =measurement,color=variable)) + 
#   ggtitle("Output (estimation)"))
#   )
# 
# p = c(p,list(
#   ggplot(data=filter(df_dataset,variable %in% c("yv","yvor"))) + 
#   geom_line(aes(x = time,y =measurement,color=variable)) + 
#   ggtitle("Output (validation)"))
#   )

p = c(p,list(
  ggplot(data=filter(df_all,variable %in% c("ye","ye_osa","ye_fr"))) + 
  geom_line(aes(x = time,y =measurement,color=variable)) + 
  ggtitle(paste0("predictions (estimation) - R2 osa = ",R2e_osa," R2 fr = ",R2e_fr)))
  )

p = c(p,list(
  ggplot(data=filter(df_all,variable %in% c("yv","yv_osa","yv_fr"))) + 
  geom_line(aes(x = time,y =measurement,color=variable)) +  
  ggtitle(paste0("predictions (validation) - R2 osa = ",R2v_osa," R2 fr = ",R2v_fr)))
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
  





