# Helon, 13/3/18
# baseado no exemplo de https://github.com/kylehamilton/deep-learning-with-r-notebooks/blob/master/notebooks/3.7-predicting-house-prices.Rmd
# Package keras version 2.1.4

rm(list=ls(all=TRUE))
cat("\014")  

library(keras)
library(ggplot2)
library(R.matlab)
source("fun_sys_id.R")
source("plot_xcorrel.R")
source("multiplot.R")

# order for the regression matrix
ny = c(1,2)
nu = c(1,2)
maxn = max(c(ny,nu))

# estudo de caso do motor 
# https://www.sciencedirect.com/science/article/pii/S127096381500070X?via%3Dihub
auxu = readMat("uval_multisine.mat")
auxy = readMat("yval_multisine.mat")
y = auxy$yval.multisine[1,]
u = auxu$uval.multisine[1,]
ndata = length(y)

y <- scale(y, center = mean(y), scale = sd(y))
u <- scale(u, center = mean(u), scale = sd(u))

# ratio for training and test
trPct = 0.5
ntr = round(trPct*ndata)

ytr = y[1:ntr]
utr = u[1:ntr]
yte = y[(ntr+1):ndata]
ute = u[(ntr+1):ndata]

# create reg matrix
list[all_data, all_targets]     = create_reg_matrix(y,u,ny,nu)
list[train_data, train_targets] = create_reg_matrix(ytr,utr,ny,nu)
list[test_data, test_targets]   = create_reg_matrix(yte,ute,ny,nu)

uxcorr = u[(maxn+1):ntr]  # u for corr tests

# # start training : k-fold CV
# k <- 4
# indices <- sample(1:nrow(train_data))
# folds <- cut(1:length(indices), breaks = k, labels = FALSE) 
# num_epochs <- 300
# all_mae_histories <- NULL
# for (i in 1:k) {
#   cat("processing fold #", i, "\n")
#   # Prepare the validation data: data from partition # k
#   val_indices <- which(folds == i, arr.ind = TRUE) 
#   val_data <- train_data[val_indices,]
#   val_targets <- train_targets[val_indices]
#   
#   # Prepare the training data: data from all other partitions
#   partial_train_data <- train_data[-val_indices,]
#   partial_train_targets <- train_targets[-val_indices]
#   
#   # Build the Keras model (already compiled)
#   model <- build_model()
#   
#   # Train the model
#   history <- model %>% fit(
#     partial_train_data, partial_train_targets,
#     validation_data = list(val_data, val_targets),
#     epochs = num_epochs, batch_size = 16, verbose = 1
#   )
#   mae_history <- history$metrics$val_mean_absolute_error
#   all_mae_histories <- rbind(all_mae_histories, mae_history)
# }
# 
# average_mae_history <- data.frame(
#   epoch = seq(1:ncol(all_mae_histories)),
#   validation_mae = apply(all_mae_histories, 2, mean)
# )
# 
# ggplot(average_mae_history, aes(x = epoch, y = validation_mae)) + geom_line() + geom_smooth()

# Get a fresh, compiled model.
model <- build_model(lr = 0.001)

# Train it on the entirety of the data.
model %>% fit(train_data, train_targets,
              epochs = 100, batch_size = 16, verbose = 1)

# result <- model %>% evaluate(test_data, test_targets) # pq deu pau?

# PREDICTIONS -------
# OSA
yh       = predict_on_batch(model, x = all_data)
yh_train = predict_on_batch(model, x = train_data)
yh_test  = predict_on_batch(model, x = test_data)

R2   = calc_R2(all_targets,yh[,]) 
R2tr = calc_R2(train_targets,yh_train[,])
R2te = calc_R2(test_targets,yh_test[,])

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

