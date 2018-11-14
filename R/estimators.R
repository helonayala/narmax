estimateNotImplemented = function (model) {
  stop(sprintf('Method "estimate" not yet implemented for model %s', class(model)))
}

#' @export
estimate = function (model, ...) UseMethod('estimate')

#' @export
estimate.default = function (model, ...) {
  stop(sprintf('Invalid class (%s) for model estimation', class(model)))
}

#' @title Estimate AR model
#' @export
estimate.ar = function (model, Y, niter = 10) {

  ny = model$ny

  Phie = genRegMatrix(model,Y)$P
  Ye   = genTarget(model,Y)

  # estimate parameters -----------------------------------------------------
  theta = MASS::ginv(Phie) %*% Ye

  model$coefficients = theta

  model$terms = colnames(Phie)

  return(model)
}

#' @title Estimate ARX model
#' @export
estimate.arx = function (model, Y, U, niter = 10) {

  ny = model$ny
  nu = model$nu

  Phie = genRegMatrix(model,Y,U)$P
  Ye   = genTarget(model,Y)

  # estimate parameters -----------------------------------------------------
  theta = MASS::ginv(Phie) %*% Ye

  model$coefficients = theta

  model$terms = colnames(Phie)

  return(model)
}


#' @title Estimate ARMA model
#' @export
estimate.arma = function (model, Y, niter = 10) {

  ny = model$ny
  ne = model$ne
  p = model$maxLag

  modelAR = ar(ny)

  Phie   = genRegMatrix(modelAR,Y)$P
  YeAR   = genTarget(modelAR,Y)
  YeARMA = genTarget(model,Y)

  # estimate parameters -----------------------------------------------------
  Th_ARMA_hat = matrix(0,ny+ne,niter)
  th_AR_hat = MASS::ginv(Phie) %*% YeAR
  th_AR_hat0 = th_AR_hat
  th_ARMA_hat0 = c(th_AR_hat0,rep(0,ne))
  ee_s1 = c(rep(0,modelAR$maxLag-1),YeAR - (Phie %*% th_AR_hat))
  dlt1 = rep(0,niter)
  dlt2 = rep(0,niter)

  for (i in 1:niter){

    Phie_ext = genRegMatrix(model,Y,ee_s1)$P

    th_ARMA_hat = MASS::ginv(Phie_ext) %*% YeARMA

    # --- stop conditions
    dlt1[i] = sqrt(sum((th_ARMA_hat - th_ARMA_hat0)^2))
    th_ARMA_hat0 = th_ARMA_hat

    # calculate error and pad zeros for the initial conditions
    ee_s = ee_s1
    ee_s1 = c(rep(0,p-1),YeARMA - (Phie_ext %*% th_ARMA_hat)[,])
    dlt2[i] = sqrt(sum((ee_s1 - ee_s)^2))

    # save estimated vectors
    Th_ARMA_hat[,i] = th_ARMA_hat
  }

  theta = th_ARMA_hat[,]

  model$coefficients = theta

  model$terms = colnames(Phie_ext)

  return(model)
}

#' @title Estimate ARMAX model
#' @export
estimate.armax = function (model, Y, U, niter = 10) {

  ny = model$ny
  nu = model$nu
  ne = model$ne
  p  = model$maxLag

  modelARX = arx(ny,nu)

  Phie = genRegMatrix(modelARX,Y,U)$P
  YeARX   = genTarget(modelARX,Y)
  YeARMAX = genTarget(model,Y)

  # estimate parameters -----------------------------------------------------
  Th_ARMAX_hat = matrix(0,ny+nu+ne,niter)
  th_ARX_hat = MASS::ginv(Phie) %*% YeARX
  th_ARX_hat0 = th_ARX_hat
  th_ARMAX_hat0 = c(th_ARX_hat0,rep(0,ne))
  ee_s1 = c(rep(0,modelARX$maxLag-1),YeARX - (Phie %*% th_ARX_hat))
  dlt1 = rep(0,niter)
  dlt2 = rep(0,niter)

  for (i in 1:niter){

    Phie_ext = genRegMatrix(model,Y,U,ee_s1)$P

    th_ARMAX_hat = MASS::ginv(Phie_ext) %*% YeARMAX

    # --- stop conditions
    dlt1[i] = sqrt(sum((th_ARMAX_hat - th_ARMAX_hat0)^2))
    th_ARMAX_hat0 = th_ARMAX_hat

    # calculate error and pad zeros for the initial conditions
    ee_s = ee_s1
    ee_s1 = c(rep(0,p-1),YeARMAX - (Phie_ext %*% th_ARMAX_hat)[,])
    dlt2[i] = sqrt(sum((ee_s1 - ee_s)^2))

    # save estimated vectors
    Th_ARMAX_hat[,i] = th_ARMAX_hat
  }

  theta = th_ARMAX_hat[,]

  model$coefficients = theta

  model$terms = colnames(Phie_ext)

  return(model)
}

#' @title Estimate NARX model
#' @export
estimate.narx = function (model, Y, U, rho = 1e-2) {

  P = genRegMatrix(model, Y, U)$P
  Target = genTarget(model, Y)

  resultNarx = frols(P, Target, rho)

  model$terms = colnames(resultNarx$Psel)
  model$coefficients = resultNarx$th

  return(model)
}

#' @title Estimate NARMA model
#' @export
estimate.narma = function (model, Y, rho_p = 1e-2, rho_n = 1.9e-6) {
  # Identify process sub-model (NAR) FROLS

  mdlNAR = nar(model$ny, model$nl)

  Pp = genRegMatrix(mdlNAR, Y)$Pp
  Target = genTarget(model, Y)

  resultNar = frols(Pp, Target, rho_p)
  E = c(rep(0, mdlNAR$maxLag - 1), Target - resultNar$W %*% resultNar$g)

  Reg = cbind(resultNar$Psel, genRegMatrix(model, Y, E = E)$Pnp)
  resultNarma = frols(Reg, Target, rho_n)

  model$terms = colnames(resultNarma$Psel)

  nth = length(resultNarma$th)
  E = c(rep(0, model$maxLag - 1), Target - resultNarma$W %*% resultNarma$g)

  iterELS = 10
  thNarmaHat = matrix(0, nrow = nth, ncol = iterELS)
  dlt = rep(0, iterELS)
  theta = NULL

  for (s in 1:iterELS) {
    P = genRegMatrix(model, Y, E=E)$P
    theta = MASS::ginv(P) %*% Target

    E1 = E
    E = c(rep(0, model$maxLag - 1), Target - (P %*% theta)[, ])

    dlt[s] = sqrt(sum((E - E1) ^ 2))
    thNarmaHat[, s] = theta
  }

  model$coefficients = theta[,]

  return(model)
}

#' @title Estimate NARMAX model
#' @export
estimate.narmax = function (model, Y, U, rho_p = 1e-2, rho_n = 1.9e-6) {
  # Identify process sub-model (NARX) FROLS
  E = rep(0, length(Y))
  Pp = genRegMatrix(model, Y, U, E)$Pp
  Target = genTarget(model, Y)

  resultNarx = frols(Pp, Target, rho_p)
  E = c(rep(0, model$maxLag - 1), Target - resultNarx$W %*% resultNarx$g)

  Reg = cbind(resultNarx$Psel, genRegMatrix(model, Y, U, E)$Pnp)
  resultNarmax = frols(Reg, Target, rho_n)

  model$terms = colnames(resultNarmax$Psel)

  nth = length(resultNarmax$th)
  E = c(rep(0, model$maxLag - 1), Target - resultNarmax$W %*% resultNarmax$g)

  iterELS = 10
  thNarmaxHat = matrix(0, nrow = nth, ncol = iterELS)
  dlt = rep(0, iterELS)
  theta = NULL

  for (s in 1:iterELS) {
    P = genRegMatrix(model, Y, U, E)$P
    theta = MASS::ginv(P) %*% Target

    E1 = E
    E = c(rep(0, model$maxLag - 1), Target - (P %*% theta)[, ])

    dlt[s] = sqrt(sum((E - E1) ^ 2))
    thNarmaxHat[, s] = theta
  }

  model$coefficients = theta[,]

  return(model)
}

#' @title Estimate ANN-ts model
#' @export
estimate.annts = function (model, Y, lr = 1e-3, epochs = 100, batch_size = 32, verbose = 1) {

  model$mdl %>% keras::compile(
    #optimizer = keras::optimizer_adam(lr=lr,amsgrad = TRUE,clipnorm=50),
    optimizer = keras::optimizer_rmsprop(lr=lr),
    loss = "mse",
    metrics = c("mae")
  )

  train_data = genRegMatrix(model,Y)
  train_targets = genTarget(model,Y)[,]

  model$mdl %>% keras::fit(train_data, train_targets, epochs = epochs, batch_size = batch_size, verbose = verbose)

#' @title Estimate caret-NARX model
#' @export
estimate.caret = function (model, Y, U, trControl = NULL,tuneGrid = NULL, tuneLength = NULL) {

  phi = data.frame(genRegMatrix(model,Y,U)$P)
  Y   = genTarget(model,Y)

  model$mdl = caret::train(phi,Y[,1],
                           method = model$method,
                           trControl = trControl,
                           tuneGrid = tuneGrid,
                           tuneLength = tuneLength,
                           verbose = TRUE,
                           metric = "Rsquared")

  # if (is.null(grid)){
  #   model$mdl = caret::train(phi,Y[,1],
  #                            method = model$method,
  #                            trControl = trControl,
  #                            tuneLength = 100,
  #                            verbose = TRUE,m)
  # } else {
  #   model$mdl = caret::train(phi,Y[,1],
  #             method = model$method,
  #             trControl = trControl,
  #             tuneGrid = tuneGrid,
  #             tuneLength = 5,
  #             verbose = TRUE)
  # }
  return(model)
}

