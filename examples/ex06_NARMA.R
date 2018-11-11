# time-series prediction example
clearWorkspace()
set.seed(42) # allows reproducibility
library(narmax)

y  = as.numeric(read.csv("../data/data_chen.csv",header=FALSE)$V1)
y = (y - min(y))/(max(y)-min(y))
n  = length(y)

# define model structure --------------------------------------------------
rho_p = 1e-5
rho_n = 1e-6

ny = 15
ne = 5
nl = 2

mdl = narma(ny,ne,nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl, y, rho_p, rho_n)
print(mdl)

# calculate predictions ---------------------------------------------------
P0 = predict(mdl, y, u, K = 0)
P1 = predict(mdl, y, u, K = 1)

# plot predictions/residuals ----------------------------------------------
print(Pe0$ploty)
print(Pe0$plote)
print(Pe1$xcorrel) # validate with correlation-based tests





mdl = annts(1:13,c(50),"sigmoid") # order from Chen paper

mdl = estimate(mdl,y2, lr = 1e-4, epochs = 100, batch_size = 32, verbose = 1)

# perform many predictions
nsteps = 1
Y = list()
g = list()
R2 = rep(0,nsteps)
for (i in 1:nsteps) {
  Y[[i]] = predict(mdl,y2,K = i) # i-step  ahead
  R2[i] = calcR2(Y[[i]]$y,Y[[i]]$yh)
  g[[i]] = ggplot(Y[[i]] %>% gather(variable,measurement,-time)) +
    geom_line(aes(x=time,y=measurement,color=variable)) +
    ggtitle(paste(i,"-steps ahead R2",R2[i]))
}

g

checkresiduals(Y[[1]]$e)
