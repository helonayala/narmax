library(narmax)

data <- read.csv('data/DATAUNIF.csv',stringsAsFactors=FALSE,header = TRUE)
u <- data$u12
y <- data$z12

mdl <- narmax(ny = 10, nu = 10, ne = 2, nl = 2)
mdl <- estimate(mdl, y, u, rho_p = 1e-3, rho_n = 1e-5)
mdl
P1 <- predict(mdl, y, u, K = 1)
P0 <- predict(mdl, y, u, K = 0)

print(P0$ploty)

print(P1$xcorrel) # validate with correlation-based tests

