
library(narmax)
library(tidyverse)
library(caret)

clearWorkspace()

# library(doParallel)
# cl <- makeCluster(parallel::detectCores())
# registerDoParallel(cl)

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
  "svmRadial",
  "svmPoly",
  "enet",
  "gamSpline",
  "glmStepAIC",
  "gaussprRadial",
  "plsRglm",
  "lm"
  )
# no good results (for this dataset and using random parameters):
# "foba"
# "xgboostlinear"
# "rvmPoly"
# "bayesglm"
# "gamboost"
# "randomGLM"
# "xgbDART"
# "glmnet"
# "kknn"
# "logreg"
# "mlpKerasDropout"
# "parRF"
# "rbf"
# "dnn"
# gbm
# bridge
# "DENFIS"
# "glm.nb"
# "neuralnet"
# qrnn
# krlsRadial
# SBC
# "avNNet",

# next we use a series of lapply functions.
# it will apply the functions to all elements in the 1st argument

# define models
mdl  = lapply(mdlOptions,caret,oy = 1:10, ou = 1:10)
# estimate models
mdl  = lapply(mdl,estimate,Y = y_est,U = u_est,trControl = trnCtrl, tuneLength = 10)
# perform predictions
Pe1  = lapply(mdl,predict,y = y_est,u = u_est, K = 1)
Pa0  = lapply(mdl,predict,y = y,    u = u,     K = 0)
# vectorize
R2e1 = as.numeric(sapply(Pe1, `[`, "R2"))
R2a0 = as.numeric(sapply(Pa0, `[`, "R2"))

# plot all R2 obtained
dfR2 = data.frame(mdlOptions,R2e1,R2a0) %>% gather(variable,value,-mdlOptions)

p = ggplot(dfR2,aes(x=mdlOptions,y=value,color=variable)) +
  geom_point()

plot(p)
