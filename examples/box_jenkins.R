
# Clear console
cat("\014")
# Clean workspace
rm(list=ls())


library(tidyverse)
library(narmax)
library(caret)

set.seed(825)

u = scan("../data/ts22.dat")
y = scan("../data/ts23.dat")

smax_inp = 1.2*max(u)
smax_out = 1.2*max(y)

y = y/smax_out
u = u/smax_inp

N = length(u) # total amount of data
Ne = 100      # total amount of data (estimation)
ind_est = 1:Ne

y_est = y[ind_est]
y_val = y[-ind_est]

u_est = u[ind_est]
u_val = u[-ind_est]

# # Step 1: Define de model
# ny = 1:3
# nu = 1:2
# mdl = caret(ny,nu,method = "svmRadialCost")
#
# model_dummy = ann(ny,nu,5,"sigmoid")
#
# phi = genRegMatrix(model_dummy,y_est,u_est)$P
#
# Y = genTarget(model_dummy,y_est)
#
# fitControl <- trainControl(
#   ## 10-fold CV
#   method = "repeatedcv",
#   number = 10,
#   ## repeated ten times
#   repeats = 10)
#
# df = data.frame(cbind(phi,Y))
#
# mdl = train(y.k. ~ ., data = df,
#             method = "svmRadial",
#             trControl = fitControl,
#             verbose = TRUE)
#
# Yse = predict(mdl, df2) # OSA

# Step 1: Define de model
# SVMPOLY
# ny = 1:3
# nu = 1:2
# mdl = caret(ny,nu,method = "svmPoly")
# trnCtrl = trainControl(number = 1,
#                        verboseIter = TRUE)
# mdl = estimate(mdl,y_est,u_est,trnCtrl)

# SVMRADIAL
# ny = 1:3
# nu = 1:2
#
# trnCtrl = trainControl(number = 5,verboseIter = TRUE) # number is how many resamples
# svmGrid = expand.grid(sigma = seq(1e-4,1e0,length.out =  10),
#                       C     = seq(1e-4,1e0,length.out =  10))
# mdl = caret(ny,nu,method = "svmRadial")
# mdl = estimate(mdl,y_est,u_est,trControl = trnCtrl,grid = svmGrid)

# KNN
# ny = 1:3
# nu = 1:2
# mdl = caret(ny,nu,method = "kknn")
# trnCtrl <- trainControl(method = "cv", number = 3, returnResamp = "all")
# mdl = estimate(mdl,y_est,u_est,trControl = trnCtrl)

# GBM
ny = 1:3
nu = 1:2
mdl = caret(ny,nu,method = "gbm")
trnCtrl <- trainControl(number = 20,
                        verboseIter = TRUE)
gbmGrid = expand.grid(interaction.depth = (1:5) * 2,
                       n.trees = (1:10)*25,
                       shrinkage = .1,
                       n.minobsinnode = 10)
gbmFit = estimate(mdl,y_est,u_est,trControl = trnCtrl,grid = gbmGrid)

print(mdl)

Yse = predict(mdl, y_est, u_est, K = 0)
Ype = predict(mdl, y_est, u_est, K = 1)
Ysv = predict(mdl, y_val, u_val, K = 0)
Ypv = predict(mdl, y_val, u_val, K = 1)

p = mdl$maxLag
time = p:N

R2se = calcR2(y_est[p:Ne],  Yse)
R2pe = calcR2(y_est[p:Ne],  Ype)
R2sv = calcR2(y_val[p:(N-Ne)],Ysv)
R2pv = calcR2(y_val[p:(N-Ne)],Ypv)

Ex = y_est[p:Ne] - Ype
Ux = u_est[p:Ne]

print(c("R2pe","R2se","R2pv","R2sv"))
print(c(R2pe,R2se,R2pv,R2sv))

dfall = data.frame(time = 1:N,
                y = y,
                u = u) %>% gather(variable, measurement, -time)

dfe = data.frame(time = p:Ne,
                 Yse = Yse,
                 Ype = Ype) %>% gather(variable, measurement, -time)

dfv = data.frame(time = (Ne+p):N,
                 Ysv = Ysv,
                 Ypv = Ypv) %>% gather(variable, measurement, -time)

df = rbind(dfall,dfe,dfv)

# predictions
p1 = ggplot(filter(df, variable %in% c("y","Yse","Ype","Ysv","Ypv"))) +
  geom_line(aes(x = time,y = measurement,color=variable)) +
  labs(title = "Predictions\n", x = "Sample", y = "Output", color = "\n")

print(p1)

g = xcorrel(Ex,Ux)

multiplot(g,cols = 1)

