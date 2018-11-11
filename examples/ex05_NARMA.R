# example from the textbook. Billings, S.A. 2013
clearWorkspace()
set.seed(42) # allows reproducibility
library(narmax)

# generate simulation data ------------------------------------------------
N = 400
e = rnorm(N, mean = 0, sd=0.04^2)
y = c(0,0,0)
for (k in 4:N) {
  y[k] = 0.5*y[k-1] +
    -0.5*y[k-2] +
    0.5*y[k-3] +
    0.01*e[k-1] +
    0.1*y[k-2]*e[k-2] +
    0.01*e[k]
}
plot(y)
lines(y)

plot(e)
lines(e)


# define model structure --------------------------------------------------
rho_p = 1e-1
rho_n = 1e-3

ny = 3
ne = 2
nl = 2

mdl = narma(ny, ne, nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl,y,rho_p,rho_n)
print(mdl)

# # calculate predictions ---------------------------------------------------
#P0 = predict(mdl, y, u, K = 0)
P1 = predict(mdl, y, u, K = 1)
#
# # plot predictions/residuals ----------------------------------------------
# print(P0$ploty)
# print(P0$plote)
# print(P1$xcorrel) # validate with correlation-based tests
#


