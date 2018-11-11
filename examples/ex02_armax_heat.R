# example from the textbook. Billings, S.A. 2013
library(narmax)
clearWorkspace()

# read input/output data --------------------------------------------------
data = matrix(scan("../data/heating_system.dat"),
               nrow=801,
               byrow=TRUE)
N = nrow(data) # number of samples
u = data[,2]
y = data[,3]
ie = 1:400
iv = 401:801
ye = y[ie]
ue = u[ie]
yv = y[iv]
uv = u[iv]

# model parameters
na = 4
nb = 4
nc = 10
mdl = armax(na,nb,nc)

# estimate parameters -----------------------------------------------------
mdl = estimate(mdl,ye,ue)
mdl

# calculate predictions ---------------------------------------------------
# estimation phase
Pe0 = predict(mdl, ye, ue, K = 0) # free-run
Pe1 = predict(mdl, ye, ue, K = 1) # one-step-ahead
# validation phase
Pv0 = predict(mdl, yv, uv, K = 0) # free-run
Pv1 = predict(mdl, yv, uv, K = 1) # one-step-ahead
# all phases
Pa0 = predict(mdl, y, u, K = 0) # free-run
Pa1 = predict(mdl, y, u, K = 1) # one-step-ahead
# print plots -------------------------------------------------------------
print(Pa0$plote)   # plot residuals
print(Pa0$ploty)   # plot predictions vs measured data
print(Pe1$xcorrel) # validate with correlation-based tests

