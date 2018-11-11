# example from the textbook. Billings, S.A. 2013
library(narmax)
clearWorkspace()

# generate simulation data ------------------------------------------------
N = 400
u = rnorm(N, mean = 0, sd=1)
e = rnorm(N, mean = 0, sd=0.04^2)
InitE = e
y = rep(0, length(u))
for (k in 3:N) {
  y[k] = 0.5*y[k-1] +
    u[k-2] +
    0.1*(u[k-2]^2) +
    0.5*e[k-1] +
    0.1*u[k-1]*e[k-2] +
    e[k]
}
e = rep(0,N)

# define model structure --------------------------------------------------
rho_p = 1e-2
rho_n = 1.9e-6

nu = 2
ny = 2
ne = 2
nl = 2

mdl = narmax(ny, nu, ne, nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl, y, u, rho_p, rho_n)
print(mdl)

# calculate predictions ---------------------------------------------------
P0 = predict(mdl, y, u, K = 0)
P1 = predict(mdl, y, u, K = 1)

# plot predictions/residuals ----------------------------------------------
print(P0$ploty)
print(P0$plote)
print(P1$xcorrel) # validate with correlation-based tests



