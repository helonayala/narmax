# Box-Jenkins furnace example
library(narmax)
clearWorkspace()

u = scan("../data/ts22.dat")
y = scan("../data/ts23.dat")

smax_inp = 1.2*max(u)
smax_out = 1.2*max(y)

y = y/smax_out
u = u/smax_inp

N  = length(u) # total amount of data
Ne = 250       # total amount of data (estimation)
ind_est = 1:Ne

ye = y[ind_est]
yv = y[-ind_est]

ue = u[ind_est]
uv = u[-ind_est]

# define model structure --------------------------------------------------
rho_p = 1e-5
rho_n = 1e-6

nu = 5
ny = 2
ne = 1
nl = 2

mdl = narmax(ny, nu, ne, nl)

# estimate the model ------------------------------------------------------
mdl = estimate(mdl, ye, ue, rho_p, rho_n)
print(mdl)

# calculate predictions ---------------------------------------------------
Pe0 = predict(mdl, ye, ue, K = 0)
Pe1 = predict(mdl, ye, ue, K = 1)
Pa0 = predict(mdl, y,  u,  K = 0)

# plot predictions/residuals ----------------------------------------------
print(Pe0$ploty)
print(Pe0$plote)
print(Pe1$xcorrel) # validate with correlation-based tests


