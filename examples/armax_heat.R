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
library(tidyverse) # data utils

# load functions
source("library_sysid.R")

data <- matrix(scan("heating_system.dat"),
               nrow=801,
               byrow=TRUE)

# allows reproducibility
set.seed(42)
#model parameters
na = 2
nb = 2
nc = 2
p  = 1 + max(na,nb,nc)

# generate input-output data ----------------------------------------------
N = nrow(data) # number of samples
u = data[,2]
y = data[,3]
ye = y[1:400]
ue = u[1:400]
yv = y[401:801]
uv = u[401:801]

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
ee_s1 = c(rep(0,p-1),Ye - (Phie %*% th_ARX_hat))
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

df_datasete = tibble(time = 1:400,
                    ye = ye,
                    ue = ue) %>%
  gather(variable, measurement, -time)

df_datasetv = tibble(time = 401:801,
                     yv = yv,
                     uv = yv) %>%
  gather(variable, measurement, -time)

df_dataset = bind_rows(df_datasete,df_datasetv)


df_prede = tibble(time = p:400,
                  ye_osa = ye_osa,
                  ye_fr  = ye_fr) %>%
  gather(variable, measurement, -time)
                  
df_predv = tibble(time = (400+p):801,
                 yv_osa = yv_osa,
                 yv_fr  = yv_fr) %>%
  gather(variable, measurement, -time)

df_pred = bind_rows(df_prede,df_predv)

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

p = ggplot(data=dplyr::filter(df_all,variable %in% c("ye","ye_osa","ye_fr","yv","yv_osa","yv_fr"))) + 
  geom_line(aes(x = time,y =measurement,color=variable))

print(p)

print(paste0("predictions (estimation) - R2 osa = ",R2e_osa," R2 fr = ",R2e_fr))
print(paste0("predictions (validation) - R2 osa = ",R2v_osa," R2 fr = ",R2v_fr))






