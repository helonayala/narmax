# example from Chen 1991 paper
library(narmax)
clearWorkspace()

y = scan("../../data/data_chen.csv")

# define model structure --------------------------------------------------
rho_p = 1e-2
rho_n = 9e-3

ny = 20
ne = 5
nl = 2
mdl = narma(ny, ne, nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl,y,rho_p,rho_n)
print(mdl)

# calculate predictions ---------------------------------------------------
P1  = predict(mdl, y, K = 1)
P2  = predict(mdl, y, K = 2)
P3  = predict(mdl, y, K = 3)
P10 = predict(mdl, y, K = 10) # takes a while to run ( we welcome parallelization :) )

# plot predictions/residuals ----------------------------------------------
print(P1$ploty)
print(P2$ploty)
print(P10$ploty)
print(P1$xcorrel) # validate with correlation-based tests


