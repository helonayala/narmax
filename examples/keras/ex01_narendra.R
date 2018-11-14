library(narmax)
library(tidyverse)

clearWorkspace()

set.seed(0)

# load data from the bouc wen benchmark
load("../../data/narendra.RData")

ne = length(ue) # amount of data - estimation
nv = length(uv) # amount of data - validation

# normalize data: mean = 0, sd = 1
y = scale(c(ye,yv))
u = scale(c(ue,uv))
ye1 = y[1:800]
yv1 = y[801:1600]
ue1 = u[1:800]
uv1 = u[801:1600]

# select model orders
oy = 1:3
ou = 1:2

# model parameters
nrn = c(128,128,128,128)
#acf = "tanh"
acf = "sigmoid"
mdl = ann(oy,ou,nrn,acf) # create model variable
#mdl = narx(3,2,2)

# estimate model parameters
mdl = estimate(mdl,ye1,ue1,lr = 1e-4, epochs = 200, batch_size = 128, verbose = 1)
#mdl = estimate(mdl,ye1,ue1,rho = 1e-3)

# PREDICITONS - estimation phase
Pe1 = predict(mdl,ye1,ue1,K = 1) # one-step-ahead
Pe0 = predict(mdl,ye1,ue1,K = 0) # free-run
Pa0 = predict(mdl,y,  u,  K = 0) # free-run

plot(Pe1$xcorrel)
plot(Pa0$ploty)


