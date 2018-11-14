
library(narmax)
# library(forecast) # for some data analysis
set.seed(123)
clearWorkspace()

y  = as.numeric(read.csv("../../data/data_chen.csv",header=FALSE)$V1)
y2 = (y - min(y))/(max(y)-min(y))
n  = length(y2)

# ggAcf(y2, lag.max = 50)
# ggAcf(diff(y2), lag.max = 50)
# ts.plot(y2)
# gglagplot(y2,lags = 30)

mdl = annts(1:18,c(75,75),"sigmoid") # order from Chen paper: 1:13
mdl = estimate(mdl,y2, lr = 1e-3, epochs = 100, batch_size = 2, verbose = 1)
# perform many predictions
P1 = predict(mdl,y2,K = 1)
P2 = predict(mdl,y2,K = 2)
P3 = predict(mdl,y2,K = 3)

plot(P1$ploty)
plot(P1$xcorrel)
plot(P3$ploty)
