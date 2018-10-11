# armax system identification
# helon - 4/9/18
# mec2015 - system identification - puc-rio

library(tidyverse)
# library(MASS)
# source("library_sysid.R")
set.seed(0) # reproducibility

# code begins

# Generate simulation data
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

# Step 1: Define de model
rho_p = 1e-2
rho_n = 1.9e-6

nu = 2
ny = 2
ne = 2
nl = 2

mdl = narmax(ny, nu, ne, nl)

# Step 2: Estimate the model
mdl = estimate(mdl, y, u, rho_p, rho_n)
print(mdl)

Ys = predict(mdl, y, u, K = 0)
Yp = predict(mdl, y, u, K = 1)

p = mdl$maxLag
time = p:N
Ep = y[p:N] - Yp
Es = y[p:N] - Ys
Up = u[p:N]

df = tibble(time,Y=y[p:N],Yp,Ys,Ep,Es) %>% gather(variable, measurement, -time)

head(df)

# predictions
p1 = ggplot(filter(df, variable %in% c("Y","Ys"))) +
  geom_line(aes(x = time,y = measurement,color=variable)) +
  labs(title = "Free-run simulation\n", x = "Sample", y = "Output", color = "\n") +
  scale_color_manual(labels = c("Measurement", "Prediction"), values = c("black", "blue")) +
  theme(legend.position="bottom")

p2 = ggplot(filter(df, variable %in% c("Y","Yp"))) +
  geom_line(aes(x = time,y = measurement,color=variable)) +
  labs(title = "One-step-ahead prediction\n", x = "Sample", y = "Output", color = "\n") +
  scale_color_manual(labels = c("Measurement", "Prediction"), values = c("black", "blue")) +
  theme(legend.position="none")

multiplot(p1,p2)

# residuals
p3 = ggplot(filter(df, variable %in% c("Es"))) +
  geom_line(aes(x = time,y = measurement)) +
  labs(title = "Free-run simulation error\n", x = "Sample", y = "Error")

p4 = ggplot(filter(df, variable %in% c("Ep"))) +
  geom_line(aes(x = time,y = measurement)) +
  labs(title = "One-step-ahead error\n", x = "Sample", y = "Error")

multiplot(p3,p4)

g = xcorrel(Ep,Up)

multiplot(g)



