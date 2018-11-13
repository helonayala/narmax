# example 06 -  sunspot data from https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/SUNSPOT/
library(narmax)
clearWorkspace()

# 1749 - 2017 monthly data
A = as.matrix(read.table("../data/sunspot.data",header = FALSE)) # 2
nlc = dim(A)

y = as.vector(t(A[,-1]))#,ncol = nlc[1]*(nlc[2]-1),nrow = 1)
y = tail(y,640) # use only last 640 samples
plot(y)
lines(y)


# define model structure --------------------------------------------------
rho_p = 1e-3
rho_n = 9e-4

ny = 3
ne = 3
nl = 2
mdl = narma(ny, ne, nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl,y,rho_p,rho_n)
print(mdl)

# calculate predictions ---------------------------------------------------
P1  = predict(mdl, y, K = 1)
P2  = predict(mdl, y, K = 2)
P5  = predict(mdl, y, K = 5)
# P10 = predict(mdl, y, K = 10)

# plot predictions/residuals ----------------------------------------------
print(P1$xcorrel) # validate with correlation-based tests
print(P5$ploty)

