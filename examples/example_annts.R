
# Clear plots
if(!is.null(dev.list())) dev.off()

# Clear console
cat("\014")
# Clean workspace
rm(list=ls())

library(narmax)
library(forecast) # for some data analysis
set.seed(123)

y  = as.numeric(read.csv("../data/data_chen.csv",header=FALSE)$V1)
y2 = (y - min(y))/(max(y)-min(y))
n  = length(y2)

ggAcf(y2, lag.max = 50)
ggAcf(diff(y2), lag.max = 50)
ts.plot(y2)

gglagplot(y2,lags = 30)

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

