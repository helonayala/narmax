# example 02 - data ftp://ftp.esat.kuleuven.be/pub/SISTA/data/thermic/heating_system.txt
library(narmax)

# read input/output data --------------------------------------------------
data <- matrix(scan("../data/heating_system.dat"),
               nrow=801,
               byrow=TRUE)
N <- nrow(data) # number of samples
u <- data[,2]
y <- data[,3]
ie <- 1:400
ye <- y[ie]
ue <- u[ie]

# model parameters
mdl <- armax(ny= 4,nu = 4,ne = 10)
mdl

# estimate parameters -----------------------------------------------------
mdl <- estimate(mdl,ye,ue)
mdl

# calculate predictions ---------------------------------------------------
# all phases
Pa0 <- predict(mdl, y, u, K = 0) # free-run
Pe1 <- predict(mdl, ye, ue, K = 1) # one-step-ahead

# print plots -------------------------------------------------------------
print(Pa0$plote)   # plot residuals
print(Pa0$ploty)   # plot predictions vs measured data
print(Pe1$xcorrel) # validate with correlation-based tests

