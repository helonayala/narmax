
library(narmax)
library(tidyverse)
library(caret)

clearWorkspace()

set.seed(825)

u = scan("../../data/ts22.dat")
y = scan("../../data/ts23.dat")

smax_inp = 1.2*max(u)
smax_out = 1.2*max(y)

y = y/smax_out
u = u/smax_inp

N = length(u) # total amount of data
Ne = 250      # total amount of data (estimation)
ind_est = 1:Ne

y_est = y[ind_est]
y_val = y[-ind_est]

u_est = u[ind_est]
u_val = u[-ind_est]

trnCtrl = trainControl(search = "random")
mdlOptions = c(
  )

# define models
model_type = "svmRadial"
# model_type = "svmPoly"
mdl  = caret(oy = 1:10, ou = 1:10,method = model_type)
# estimate models
mdl  = estimate(mdl,Y = y_est,U = u_est,trControl = trnCtrl, tuneLength = 100)
# perform predictions
Pe1  = predict(mdl,y = y_est,u = u_est, K = 1)
Pa0  = predict(mdl,y = y,    u = u,     K = 0)

plot(Pe1$xcorrel)
plot(Pa0$ploty)
