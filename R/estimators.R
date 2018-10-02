estimateNotImplemented = function (model) {
  stop(sprintf('Method "estimate" not yet implemented for model %s', class(model)))
}

#' @export
estimate = function (model, ...) UseMethod('estimate')

#' @export
estimate.default = function (model, ...) {
  stop(sprintf('Invalid class (%s) for model estimation', class(model)))
}

#' @title Estimate ARX model
#' @export
estimate.arx = function (model, ...) {
  estimateNotImplemented(model)
}

#' @title Estimate ARMAX model
#' @export
estimate.armax = function (model, Y, U, niter = 10) {

  ny = model$ny
  nu = model$nu
  ne = model$ne
  p = model$maxLag

  modelARX = arx(ny,nu)

  Phie = genRegMatrix(modelARX,Y,U)$P
  Ye   = genTarget(modelARX,Y)

  # estimate parameters -----------------------------------------------------
  Th_ARMAX_hat = matrix(0,ny+nu+ne,niter)
  th_ARX_hat = ginv(Phie) %*% Ye
  th_ARX_hat0 = th_ARX_hat
  th_ARMAX_hat0 = c(th_ARX_hat0,rep(0,ne))
  ee_s1 = c(rep(0,p-1),Ye - (Phie %*% th_ARX_hat))
  dlt1 = rep(0,niter)
  dlt2 = rep(0,niter)

  for (i in 1:niter){

    Phie_ext = genRegMatrix(model,Y,U,ee_s1)$P

    th_ARMAX_hat = ginv(Phie_ext) %*% Ye

    # --- stop conditions
    dlt1[i] = sqrt(sum((th_ARMAX_hat - th_ARMAX_hat0)^2))
    th_ARMAX_hat0 = th_ARMAX_hat

    # calculate error and pad zeros for the initial conditions
    ee_s = ee_s1
    ee_s1 = c(rep(0,p-1),Ye - (Phie_ext %*% th_ARMAX_hat)[,])
    dlt2[i] = sqrt(sum((ee_s1 - ee_s)^2))

    # save estimated vectors
    Th_ARMAX_hat[,i] = th_ARMAX_hat
  }

  theta = th_ARMAX_hat[,]

  model$coefficients = theta

  print(theta)

  return(model)

}

#' @title Estimate NARX model
#' @export
estimate.narx = function (model, ...) {
  estimateNotImplemented(model)
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

  print('---------------------')
  print(' NARX TERM SELECTION ')
  print('---------------------')
  print(resultNarx$ERR * 100)
  print('SUM ERR:')
  print(sum(resultNarx$ERR * 100))
  print('SELECTED TERMS:')
  print(colnames(resultNarx$Psel))

  Reg = cbind(resultNarx$Psel, genRegMatrix(model, Y, U, E)$Pnp)
  resultNarmax = frols(Reg, Target, rho_n)

  model$terms = colnames(resultNarmax$Psel)

  print(' ')
  print('-----------------------')
  print(' NARMAX TERM SELECTION ')
  print('-----------------------')
  print(resultNarmax$ERR * 100)
  print('SUM ERR:')
  print(sum(resultNarmax$ERR * 100))
  print('SELECTED TERMS:')
  print(colnames(resultNarmax$Psel))

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

  print(theta)

  return(model)
}
