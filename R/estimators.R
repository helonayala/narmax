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
estimate.armax = function (model, ...) {
  estimateNotImplemented(model)
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
