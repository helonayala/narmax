sizeGuard = function (Y, U, E = NULL) {
  if (length(Y) != length(U)) stop('Input-Output vector must have the same size')
  if (!is.null(E)) {
    if (length(Y) != length(E)) stop('Error-Output vector must have the same size')
  }
}

#' @export
generateRM = function (model, ...) {
  UseMethod('generateRM', model)
}

#' @export
generateRM.default = function (model, ...) {
  stop(sprintf('Unknown class %s for regression matrix generation', class(model)))
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for an ARX Model
#' @param model ARX Model
#' @param Y The target vector
#' @param U The input vector
#' @return Regression Matrix
#' @export
generateRM.arxModel = function (model, Y, U, E = NULL) {
  na = model$ny
  nb = model$nu
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)

  Phi = matrix(0, nrow = N - p + 1, ncol = na + nb)
  colPhi = NULL
  for(i in 1:na){
    print(i)
    print(Phi[, i])
    print(-Y[(p - i):(N - i)])
    Phi[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  for(i in 1:nb){
    Phi[, na + i] = U[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("u(k-",i,")"))
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(Phi) = colPhi
  rownames(Phi) = rowPhi

  return(Phi)
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for an ARMAX Model
#' @param model ARMAX Model
#' @param Y The target vector
#' @param U The input vector
#' @param E The error vector (can be NULL)
#' @return Regression Matrix
#' @export
generateRM.armaxModel = function (model, Y, U, E = NULL) {
  na = model$ny
  nb = model$nu
  nc = model$ne
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)

  if (is.null(E)) nc = 0

  Phi = matrix(0, nrow = N - p + 1, ncol = na + nb + nc)

  colPhi = NULL
  for(i in 1:na){
    Phi[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  for(i in 1:nb){
    Phi[, na + i] = U[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("u(k-",i,")"))
  }

  if (nc > 0) {
    for(i in 1:nc){
      Phi[, na + nb + i] = E[(p - i):(N - i)]
      colPhi = c(colPhi, paste0("e(k-",i,")"))
    }
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(Phi) = colPhi
  rownames(Phi) = rowPhi

  return(Phi)
}
