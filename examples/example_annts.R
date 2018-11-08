
# Clear console
cat("\014")
# Clean workspace
rm(list=ls())

library(narmax)

y  = as.numeric(read.csv("../data/data_chen.csv",header=FALSE)$V1)
y2 = scale(y)[,]
n  = length(y2)

mdl = annts(1:13,c(30,30,30),"sigmoid") # order from Chen paper

mdl = estimate(mdl,y2, lr = 1e-3, epochs = 20, batch_size = 1, verbose = 1)

# perform many predictions
nsteps = 5
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





