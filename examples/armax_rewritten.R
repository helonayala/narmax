# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

clearWorkspace()

# load libraries
library(ggplot2) # fancy plots
library(signal)  # filter for input signal
library(MASS)    # use ginv
library(tidyverse) # data utils
library(narmax)

# load functions
source("../organize/library_sysid.R")

# allows reproducibility
set.seed(42)
#model parameters
na = 4
nb = 3
nc = 2

mdl = armax(na,nb,nc)

# generate input-output data ----------------------------------------------
sd_noise = 1e-2
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


mdl = estimate(mdl,ye,ue)

a=2

#
#
# # calculate predictions ---------------------------------------------------
# ye_osa = calcOSA_ARMAX(ye,ue,na,nb,nc,p,th_ARMAX_hat)
# ye_fr  = calcFR_ARX(ye,ue,na,nb,p,th_ARX_hat)
# yv_osa = calcOSA_ARMAX(yv,uv,na,nb,nc,p,th_ARMAX_hat)
# yv_fr  = calcFR_ARX(yv,uv,na,nb,p,th_ARX_hat)
#
# # prediction performance
# R2e_osa = calcR2(Ye,ye_osa)
# R2e_fr  = calcR2(Ye,ye_fr)
# R2v_osa = calcR2(Yv,yv_osa)
# R2v_fr  = calcR2(Yv,yv_fr)
#
# # create tibbles for plotting ---------------------------------------------
#
# df_dataset = tibble(time = 1:N,
#                     ye = ye,
#                     yeor = yeor,
#                     ue = ue,
#                     yv = yv,
#                     yvor = yvor,
#                     uv = uv) %>%
#   gather(variable, measurement, -time)
#
# df_pred = tibble(time = p:N,
#                  ye_osa = ye_osa,
#                  ye_fr  = ye_fr,
#                  yv_osa = yv_osa,
#                  yv_fr  = yv_fr) %>%
#   gather(variable, measurement, -time)
#
# df_error = tibble(time = p:N,
#                   ee_osa = (Ye - ye_osa)[,],
#                   ee_fr  = (Ye - ye_fr)[,],
#                   ev_osa = (Yv - yv_osa)[,],
#                   ev_fr  = (Yv - yv_fr)[,]) %>%
#   gather(variable, measurement, -time)
#
# df_all = bind_rows(df_dataset,df_pred)
#
#
# # plots -------------------------------------------------------------------
# p = list()
#
# # p = c(p,list(
# #   ggplot(data=filter(df_dataset,variable %in% c("ue","uv"))) +
# #     geom_line(aes(x = time,y =measurement,color=variable)) +
# #     ggtitle("Input"))
# #   )
#
# # p = c(p,list(
# #   ggplot(data=filter(df_dataset,variable %in% c("ye","yeor"))) +
# #   geom_line(aes(x = time,y =measurement,color=variable)) +
# #   ggtitle("Output (estimation)"))
# #   )
# #
# # p = c(p,list(
# #   ggplot(data=filter(df_dataset,variable %in% c("yv","yvor"))) +
# #   geom_line(aes(x = time,y =measurement,color=variable)) +
# #   ggtitle("Output (validation)"))
# #   )
#
# p = c(p,list(
#   ggplot(data=dplyr::filter(df_all,variable %in% c("ye","ye_osa","ye_fr"))) +
#     geom_line(aes(x = time,y =measurement,color=variable)) +
#     ggtitle(paste0("predictions (estimation) - R2 osa = ",R2e_osa," R2 fr = ",R2e_fr)))
# )
#
# p = c(p,list(
#   ggplot(data=dplyr::filter(df_all,variable %in% c("yv","yv_osa","yv_fr"))) +
#     geom_line(aes(x = time,y =measurement,color=variable)) +
#     ggtitle(paste0("predictions (validation) - R2 osa = ",R2v_osa," R2 fr = ",R2v_fr)))
# )
#
# p = c(p,list(
#   ggplot(data=dplyr::filter(df_error,variable %in% c("ee_osa","ee_fr"))) +
#     geom_line(aes(x = time,y =measurement,color=variable)) +
#     ggtitle("prediction residuals (estimation)"))
# )
#
# p = c(p,list(
#   ggplot(data=dplyr::filter(df_error,variable %in% c("ev_osa","ev_fr"))) +
#     geom_line(aes(x = time,y =measurement,color=variable)) +
#     ggtitle("prediction residuals (validation)"))
# )
#
# print(p)
#





