# example 01 - simulated armax system
library(narmax)
clearWorkspace()
set.seed(42) # allows reproducibility

# define model parameters -------------------------------------------------
na = 4
nb = 3
nc = 2
mdl = armax(na,nb,nc)

# generate input-output data ----------------------------------------------
sd_noise = 1e-2
N = 1000 # number of samples
cutoff = 0.1 # normalized frequency cutoff
th = c(0.3,0.5,0.1,0.4,0.2,0.32,0.1)
# create input signal
ue = multisine(N,cutoff)
uv = randnoise(N,cutoff)

ee = rnorm(N,mean=0,sd=sd_noise)
ev = rnorm(N,mean=0,sd=sd_noise)
M_spec(ue)
# M_spec(uv,'uv')

ye = array(0,N)
yv = array(0,N)

for (k in 5:N){
  phie = c(-ye[k-1],-ye[k-2],-ye[k-3],ue[k-1],ue[k-2],ee[k-1],ee[k-2])
  phiv = c(-yv[k-1],-yv[k-2],-yv[k-3],uv[k-1],uv[k-2],ev[k-1],ev[k-2])
  ye[k] = phie %*% th
  yv[k] = phiv %*% th
}
yeor = ye
ye = ye + ee
yvor = yv
yv = yv + ev

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


