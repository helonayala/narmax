clearWorkspace()
set.seed(42) # allows reproducibility
library(narmax)

# define model parameters -------------------------------------------------
na = 4
mdl = ar(na)

# generate input-output data ----------------------------------------------
N = 1000 # number of samples
cutoff = 0.1 # normalized frequency cutoff
th = c(0.3,0.5,0.1,0.4)

y = rnorm(N,mean=0,sd=1e-2)

for (k in 5:N){
  phi = c(-y[k-1],-y[k-2],-y[k-3],-y[k-4])
  y[k] = phi %*% th
}
plot(y)

# estimate model parameters -----------------------------------------------
mdl = estimate(mdl, ye, ue)

# make model predictions --------------------------------------------------
Pe0 = predict(mdl, ye, ue, K = 0) # free-run
Pe1 = predict(mdl, ye, ue, K = 1) # one-step-ahead
Pv0 = predict(mdl, yv, uv, K = 0) # free-run
Pv1 = predict(mdl, yv, uv, K = 1) # one-step-ahead

# output plots ------------------------------------------------------------
print(Pv0$plote)   # plot residuals
print(Pv0$ploty)   # plot predictions vs measured data
print(Pe1$xcorrel) # validate with correlation-based tests

