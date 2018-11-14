sizeGuard = function (Y, U, E = NULL) {
  if (!is.null(U) && length(Y) != length(U)) stop('Input-Output vector must have the same size')
  if (!is.null(E)) {
    if (length(Y) != length(E)) stop('Error-Output vector must have the same size')
  }
}

#' @export
genRegMatrix = function (model, ...) UseMethod('genRegMatrix', model)

#' @export
genRegMatrix.default = function (model, ...) {
  stop(sprintf('Unknown class %s for regression matrix generation', class(model)))
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for an AR Model
#' @param model AR Model
#' @param Y The target vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.ar = function (model, Y, U = NULL, E = NULL) {
  na = model$ny
  p = model$maxLag
  N = length(Y)

  U = NULL
  E = NULL

  obj = list()

  obj$P = matrix(0, nrow = N - p + 1, ncol = na)
  colPhi = NULL
  for(i in 1:na){
    obj$P[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  obj$Pp = obj$P
  return(obj)
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for an ARX Model
#' @param model ARX Model
#' @param Y The target vector
#' @param U The input vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.arx = function (model, Y, U, E = NULL) {
  na = model$ny
  nb = model$nu
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)
  obj = list()

  obj$P = matrix(0, nrow = N - p + 1, ncol = na + nb)
  colPhi = NULL
  for(i in 1:na){
    obj$P[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  for(i in 1:nb){
    obj$P[, na + i] = U[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("u(k-",i,")"))
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  obj$Pp = obj$P
  return(obj)
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for an ARMAX Model
#' @param model ARMAX Model
#' @param Y The target vector
#' @param U The input vector
#' @param E The error vector (can be NULL)
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.arma = function (model, Y, U = NULL, E) {
  na = model$ny
  nc = model$ne
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)
  obj = list()

  if (is.null(E)) nc = 0

  obj$P = matrix(0, nrow = N - p + 1, ncol = na + nc)

  colPhi = NULL
  for(i in 1:na) {
    obj$P[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  if (nc > 0) {
    for(i in 1:nc) {
      obj$P[, na + i] = E[(p - i):(N - i)]
      colPhi = c(colPhi, paste0("e(k-",i,")"))
    }
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }
  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)

  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))

  obj$Pnp = matrix(
    obj$P[, errIndexes],
    ncol = sum(errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))

  return(obj)
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for an ARMAX Model
#' @param model ARMAX Model
#' @param Y The target vector
#' @param U The input vector
#' @param E The error vector (can be NULL)
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.armax = function (model, Y, U, E) {
  na = model$ny
  nb = model$nu
  nc = model$ne
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y,U, E)
  obj = list()

  if (is.null(E)) nc = 0

  obj$P = matrix(0, nrow = N - p + 1, ncol = na + nb + nc)

  colPhi = NULL
  for(i in 1:na) {
    obj$P[, i] = -Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("-y(k-",i,")"))
  }

  for(i in 1:nb) {
    obj$P[, na + i] = U[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("u(k-",i,")"))
  }

  if (nc > 0) {
    for(i in 1:nc) {
      obj$P[, na + nb + i] = E[(p - i):(N - i)]
      colPhi = c(colPhi, paste0("e(k-",i,")"))
    }
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(obj$P) = colPhi
  rownames(obj$P) = rowPhi

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }
  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)

  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))

  obj$Pnp = matrix(
    obj$P[, errIndexes],
    ncol = sum(errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))

  return(obj)
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for a NARX Model
#' @param model NAR Model
#' @param Y The target vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.nar = function (model, Y, U = NULL, E = NULL) {
  ny = model$ny
  nl = model$nl
  n = ny
  p = model$maxLag
  N = length(Y)

  obj = list()

  auxExp = list()
  candList = list()

  for (i in 1:nl) {
    eval(parse(text = paste0('auxExp$x', i, '= 1:n')))
    cand = expand.grid(auxExp)
    candList[[i]] = unique(t(apply(cand, 1, sort)))
  }

  P0 = genRegMatrix(ar(ny), Y)$P
  P0[, 1:ny] = -P0[, 1:ny]
  colnames(P0) = c(paste0("y(k-", 1:ny, ")"))

  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Constant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)

  if (nl >= 2) {
    for (i in 2:nl) {
      ncand = nrow(candList[[i]])
      for (j in 1:ncand) {
        Pcand_a = subsetMatrix(P0, NULL, candList[[i]][j, ]) # P0[, candList[[i]][j, ]]
        names = colnames(Pcand_a)
        Pcand_b = matrix(apply(Pcand_a, 1, prod), ncol = 1)
        colnames(Pcand_b) = stringr::str_c(names, collapse = '')
        obj$P = cbind(obj$P, Pcand_b)
      }
    }
  }

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  obj$Pp = obj$P
  obj$Pnp = NULL
  return(obj)
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for a NARX Model
#' @param model NARX Model
#' @param Y The target vector
#' @param U The input vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.narx = function (model, Y, U, E = NULL) {
  ny = model$ny
  nu = model$nu
  nl = model$nl
  n = ny + nu
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)
  obj = list()

  auxExp = list()
  candList = list()

  for (i in 1:nl) {
    eval(parse(text = paste0('auxExp$x', i, '= 1:n')))
    cand = expand.grid(auxExp)
    candList[[i]] = unique(t(apply(cand, 1, sort)))
  }

  P0 = genRegMatrix(arx(ny, nu), Y, U, E)$P
  P0[, 1:ny] = -P0[, 1:ny]
  colnames(P0) = c(paste0("y(k-", 1:ny, ")"), paste0("u(k-", 1:nu, ")"))

  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Contant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)

  if (nl >= 2) {
    for (i in 2:nl) {
      ncand = nrow(candList[[i]])
      for (j in 1:ncand) {
        Pcand_a = subsetMatrix(P0, NULL, candList[[i]][j, ]) # P0[, candList[[i]][j, ]]
        names = colnames(Pcand_a)
        Pcand_b = matrix(apply(Pcand_a, 1, prod), ncol = 1)
        colnames(Pcand_b) = stringr::str_c(names, collapse = '')
        obj$P = cbind(obj$P, Pcand_b)
      }
    }
  }

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  obj$Pp = obj$P
  obj$Pnp = NULL
  return(obj)
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for a NARMA Model
#' @param model NARMA Model
#' @param Y The target vector
#' @param E The error vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.narma = function (model, Y, U = NULL, E) {
  ny = model$ny
  ne = model$ne
  nl = model$nl
  n = ny + ne
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y,U,E)
  obj = list()

  auxExp = list()
  candList = list()

  for (i in 1:nl) {
    eval(parse(text = paste0('auxExp$x', i, '= 1:n')))
    cand = expand.grid(auxExp)
    candList[[i]] = unique(t(apply(cand, 1, sort)))
  }

  P0 = genRegMatrix(arma(ny, ne), Y, U, E)$P
  P0[, 1:ny] = -P0[, 1:ny]
  colnames(P0) = c(paste0('y(k-', 1:ny, ')'), paste0('e(k-', 1:ne, ')'))

  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Contant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)

  if (nl >= 2) {
    for (i in 2:nl) {
      ncand = nrow(candList[[i]])
      for (j in 1:ncand) {
        Pcand_a = subsetMatrix(P0, NULL, candList[[i]][j, ]) # P0[, candList[[i]][j, ]]
        names = colnames(Pcand_a)
        Pcand_b = matrix(apply(Pcand_a, 1, prod), ncol = 1)
        colnames(Pcand_b) = stringr::str_c(names, collapse = '')
        obj$P = cbind(obj$P, Pcand_b)
      }
    }
  }

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)
  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))

  obj$Pnp = matrix(
    obj$P[, errIndexes],
    ncol = sum(errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))

  return(obj)
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for a NARMAX Model
#' @param model NARMAX Model
#' @param Y The target vector
#' @param U The input vector
#' @param E The error vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.narmax = function (model, Y, U, E) {
  ny = model$ny
  nu = model$nu
  ne = model$ne
  nl = model$nl
  n = ny + nu + ne
  p = model$maxLag
  N = length(Y)

  sizeGuard(Y, U, E)
  obj = list()

  auxExp = list()
  candList = list()

  for (i in 1:nl) {
    eval(parse(text = paste0('auxExp$x', i, '= 1:n')))
    cand = expand.grid(auxExp)
    candList[[i]] = unique(t(apply(cand, 1, sort)))
  }

  P0 = genRegMatrix(armax(ny, nu, ne), Y, U, E)$P
  P0[, 1:ny] = -P0[, 1:ny]
  colnames(P0) = c(paste0('y(k-', 1:ny, ')'), paste0('u(k-', 1:nu, ')'), paste0('e(k-', 1:ne, ')'))

  NP0 = nrow(P0)
  obj$P = NULL
  obj$P = cbind(rep(1, NP0), obj$P) # Contant
  colnames(obj$P) = "constant"
  obj$P = cbind(obj$P, P0)

  if (nl >= 2) {
    for (i in 2:nl) {
      ncand = nrow(candList[[i]])
      for (j in 1:ncand) {
        Pcand_a = subsetMatrix(P0, NULL, candList[[i]][j, ]) # P0[, candList[[i]][j, ]]
        names = colnames(Pcand_a)
        Pcand_b = matrix(apply(Pcand_a, 1, prod), ncol = 1)
        colnames(Pcand_b) = stringr::str_c(names, collapse = '')
        obj$P = cbind(obj$P, Pcand_b)
      }
    }
  }

  if (!is.null(model$terms)) {
    obj$P = subsetMatrix(obj$P, NULL, model$terms)
  }

  errIndexes = grepl('e(', colnames(obj$P), fixed = TRUE)
  obj$Pp = matrix(
    obj$P[, !errIndexes],
    ncol = sum(!errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[!errIndexes]))

  obj$Pnp = matrix(
    obj$P[, errIndexes],
    ncol = sum(errIndexes),
    dimnames = list(rownames(obj$P), colnames(obj$P)[errIndexes]))

  return(obj)
}


#' @title Generates a regression matrix
#' @description Generates a regression matrix for a ANN-NARX Model
#' @param model ANN Model
#' @param Y The output data vector
#' @param U The input data vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.ann = function (model, Y, U, E = NULL) {
  oy = model$oy
  ou = model$ou

  ny = max(oy)
  nu = max(ou)

  sizeGuard(Y, U, E)

  P0 = genRegMatrix(arx(ny, nu), Y, U, E)$P
  P0[, 1:ny] = -P0[, 1:ny]
  colnames(P0) = c(paste0("y(k-", 1:ny, ")"), paste0("u(k-", 1:nu, ")"))

  finalReg = c(paste0("y(k-", oy, ")"), paste0("u(k-", ou, ")"))

  P = P0[,is.element(colnames(P0),finalReg)]

  obj = list()
  obj$P = P
  obj$Pp = obj$P
  obj$Pnp = NULL
  return(obj)
}

#' @title Generates a regression matrix
#' @description Generates a regression matrix for a caret-based Model
#' @param model any caret created model
#' @param Y The output data vector
#' @param U The input data vector
#' @return Object containing:
#' \describe{
#'  \item{P}{Regression matrix with all terms}
#'  \item{Pp}{Regression matrix with only process terms}
#'  \item{Pnp}{Regression matrix without process terms}
#' }
#' @export
genRegMatrix.caret = function (model, Y, U, E = NULL) {
  return(genRegMatrix.ann(model,Y,U)) # the regression matrix is the same as in ANNs
}
  
#' @title Generates a regression matrix
#' @description Generates a regression matrix for a ANN-ts Model
#' @param model ANN Model
#' @param Y The output data vector
#' @return Regression matrix:
#' \describe{
#'  \item{Regression matrix with all terms}: y(t-1), ... contained in model information
#' }
#' @export
genRegMatrix.annts = function (model, Y, U = NULL, E = NULL) {
  oy = model$oy

  ny = max(oy)

  p = model$maxLag
  N = length(Y)

  P = matrix(0, nrow = N - p + 1, ncol = ny)
  colPhi = NULL
  for(i in 1:ny){
    P[, i] = Y[(p - i):(N - i)]
    colPhi = c(colPhi, paste0("y(k-",i,")"))
  }

  rowPhi = paste0(rep("k=", N - p + 1), p:N)

  colnames(P) = colPhi
  rownames(P) = rowPhi

  finalReg = paste0("y(k-", oy, ")")

  P = P[,is.element(colnames(P),finalReg)]

  return(P)
}

#' @export
genTarget = function (model, ...) UseMethod('genTarget')

#' @title Generate target vector
#' @description Generate target vector based on maximum lag. Works with any model
#' @param model Any model containing a $maxLag property
#' @param Y Original target vector
#' @export
genTarget.default = function (model, Y) {
  N = length(Y)
  p = model$maxLag
  target = matrix(Y[p:N], ncol = 1)
  rownames(target) = paste0(rep("k=", N - p + 1), p:N)
  colnames(target) = "y(k)"
  return(target)
}

