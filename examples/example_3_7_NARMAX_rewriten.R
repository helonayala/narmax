# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

clearWorkspace()

# library(tidyverse)
# library(MASS)
# source("library_sysid.R")
set.seed(0) # reproducibility

# code begins

# Generate simulation data
N = 400
u = rnorm(N, mean = 0, sd=1)
e = rnorm(N, mean = 0, sd=0.04^2)
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

# Step 1: Define de model
rho_p = 1e-2
rho_n = 1.9e-6

nu = 2
ny = 2
ne = 2
nl = 2

model = narmax(ny, nu, ne, nl)

# Step 2: Estimate the model
model = estimate(model, y, u, rho_p, rho_n)
print(model)
