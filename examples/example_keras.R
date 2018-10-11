library(narmax)
library(tidyverse)

# load data from the bouc wen benchmark
load("../data/boucwen.RData")

# select multisine data for estimation
ue = umultisine
ye = ymultisine
uv = usinesweep
yv = ysinesweep

ne = length(ue) # amount of data - estimation
nv = length(uv) # amount of data - validation

# select model orders
oy = 1:3
ou = 1:2

# model parameters
lr = 1e-3 # learning rate
nrn = c(6,6,6)
acf = "tanh"
mdl = ann(oy,ou,nrn,acf) # create model variable

# normalize data
ye1 = scale(ye, center = mean(ye), scale = sd(ye))
ue1 = scale(ue, center = mean(ue), scale = sd(ue))
uxcorr = ue1[mdl$p:ne]  # u for corr tests

# estimate model parameters
mdl = estimate(mdl,ye1,ue1,lr)

# Train it on the entirety of the data.
model %>% fit(train_data, train_targets,
              epochs = 100, batch_size = 16, verbose = 1)











build_model <- function(lr) {
  model <- keras_model_sequential() %>%
    layer_dense(units = 8, activation = "tanh",
                input_shape = dim(train_data)[[2]]) %>%
    layer_dense(units = 8, activation = "tanh") %>%
    layer_dense(units = 8, activation = "tanh") %>%
    # layer_dense(units = 8, activation = "tanh") %>%
    # layer_dense(units = 8, activation = "tanh") %>%
    # layer_dense(units = 20, activation = "tanh") %>%
    # layer_dense(units = 20, activation = "tanh") %>%
    # layer_dense(units = 20, activation = "tanh") %>%
    layer_dense(units = 1)
  # layer_dense(units = 8, activation = "tanh",
  #             input_shape = dim(train_data)[[2]]) %>%
  #   layer_dense(units = 8, activation = "tanh") %>%
  #   layer_dense(units = 8, activation = "tanh") %>%
  #   layer_dense(units = 8, activation = "tanh") %>%
  #   layer_dense(units = 1)

  model %>% compile(
    #optimizer = optimizer_rmsprop(lr=lr),
    optimizer = optimizer_adam(lr=lr,amsgrad = TRUE,clipnorm=50),
    loss = "mse",
    metrics = c("mae")
  )
}

















# result <- model %>% evaluate(test_data, test_targets) # pq deu pau?

# PREDICTIONS -------
# OSA
yh       = predict_on_batch(model, x = all_data)
yh_train = predict_on_batch(model, x = train_data)
yh_test  = predict_on_batch(model, x = test_data)

R2   = calcR2(all_targets,yh[,])
R2tr = calcR2(train_targets,yh_train[,])
R2te = calcR2(test_targets,yh_test[,])

# FREE RUN
yh_fr       = predictFreeRun(u,y,ny,nu,model)
yh_train_fr = predictFreeRun(utr,ytr,ny,nu,model)
yh_test_fr  = predictFreeRun(ute,yte,ny,nu,model)

R2_fr   = calc_R2(all_targets,yh_fr)
R2tr_fr = calc_R2(train_targets,yh_train_fr)
R2te_fr = calc_R2(test_targets,yh_test_fr)

print(paste("R2",R2,"R2tr",R2tr,"R2te",R2te))
print(paste("R2_fr",R2_fr,"R2tr_fr",R2tr_fr,"R2te_fr",R2te_fr))
# PLOTS

plot_xcorrel(train_targets - yh_train[,], uxcorr)

ndata = length(all_targets)

df = data.frame(
  t = rep(1:ndata,3),
  y = c(all_targets,yh_fr,yh[,]),
  Type = c(rep('Measured data',ndata),
           rep('Free-run simulation',ndata),
           rep('One step ahead',ndata))
)

ggplot(df,aes(x=t,y=y,color=Type)) + geom_line() +xlab('Sample')+ylab('Output (normalized)')

# grafico de Y vs Yh
df = data.frame(
  y = all_targets,
  yh_fr = yh_fr,
  yh = yh[,]
)

p1 = ggplot(df,aes(x=y,y=yh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="blue", size=1.5)+
  xlab('Measured')+ylab('Predicted (one-step-ahead)')+ggtitle(paste("R2 =",R2))
p2 = ggplot(df,aes(x=y,y=yh_fr)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="blue", size=1.5)+
  xlab('Measured')+ylab('Predicted (free-run)')+ggtitle(paste("R2 =",R2_fr))
multiplot(p1, p2, cols = 2)

