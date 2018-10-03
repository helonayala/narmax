# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

clearWorkspace()

# load libraries
# library(ggplot2) # fancy plots
# library(signal)  # filter for input signal
# library(MASS)    # use ginv
# library(tidyverse) # data utils
# library(narmax)

# load functions
#source("../organize/library_sysid.R")

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

Yp = predict(mdl, ye, ue, K = 1)
Ys = predict(mdl, ye, ue, K = 0)

p = mdl$maxLag
df = data.frame(time = time[p:N], y = y[p:N], yp = Yp, ys = Ys)

p = ggplot(data = df, aes(x = time)) +
  geom_line(aes(y = y)) +
  geom_line(aes(y = yp), color = "red") +
  geom_line(aes(y = ys), color = "blue")

plot(p)
