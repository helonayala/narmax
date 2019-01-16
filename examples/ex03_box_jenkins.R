# Box-Jenkins furnace example
library(narmax)

u <- scan("../data/ts22.dat")
y <- scan("../data/ts23.dat")

ne <- 250
ue <- u[1:ne]
ye <- y[1:ne]

# define model structure --------------------------------------------------
mdl1 <- narmax(ny = 1, nu = 1, ne = 1, nl = 2)
mdl2 <- narmax(ny = 5, nu = 5, ne = 1, nl = 2)

# estimate the model ------------------------------------------------------
mdl1 <- estimate(mdl1, ye, ue, rho_p = 1e-5, rho_n = 1e-6)
mdl2 <- estimate(mdl2, ye, ue, rho_p = 1e-5, rho_n = 1e-6)
print(mdl1)
print(mdl2)

# calculate predictions ---------------------------------------------------
# mdl1
Pe1_1 <- predict(mdl1, ye, ue, K = 1)
Pa0_1 <- predict(mdl1, y, u, K = 0)

# mdl2
Pe1_2 <- predict(mdl2, ye, ue, K = 1)
Pa0_2 <- predict(mdl2, y, u, K = 0)

# plot predictions/residuals ----------------------------------------------
print(Pa0_1$ploty)
print(Pa0_2$ploty)

print(Pe1_1$xcorrel) # validate with correlation-based tests
print(Pe1_2$xcorrel) # validate with correlation-based tests

