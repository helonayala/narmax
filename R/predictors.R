#' @export
predict = function (model, ...) UseMethod('predict')

#' @export
predict.arx = function (model, ...) {
  cat('Running arx prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.armax = function (model, y, u, K) {
  cat('Running armax prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.narx = function (model, ...) {
  cat('Running narx prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.narmax = function (model, y, u, K = 1, ...) {
  cat('Running narmax prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.default = function (model, y, u, K, ...) {
  if (K < 0) stop('K must be greater or equal to zero')
  method = switch(
    as.character(K),
    "1" = oneStepAhead,
    "0" = freeRun,
    kStepAhead
  )
  return(method(model, y, u, K))
}

oneStepAhead = function (model, y, u, ...) {
  theta = as.matrix(model$coefficients)
  p = model$maxLag

  ySlice = y[1:(p - 1)]
  uSlice = u[1:(p - 1)]
  eSlice = e[1:(p - 1)]

  N = length(y)

  for (k in p:N) {
    auxY = c(y[(k - p):(k - 1)], 0)
    auxU = c(uSlice[(k - p):(k - 1)], 0)
    auxE = c(eSlice[(k - p):(k - 1)], 0)
    phiK = genRegMatrix(model, auxY, auxU, auxE)$P
    ySlice[k] = (phiK %*% theta)[1]
    uSlice[k] = u[k]
    eSlice[k] = y[k] - ySlice[k]
  }

  return(ySlice[p:N])
}

freeRun = function (model, y, u, K, ...) {
  theta = as.matrix(model$coefficients)
  p = model$maxLag

  ySlice = y[1:(p - 1)]
  uSlice = u[1:(p - 1)]
  eSlice = e[1:(p - 1)]

  N = length(y)

  for (k in p:N) {
    auxY = c(ySlice[(k - p):(k - 1)], 0)
    auxU = c(uSlice[(k - p):(k - 1)], 0)
    auxE = c(eSlice[(k - p):(k - 1)], 0)
    phiK = genRegMatrix(model, auxY, auxU, auxE)$P
    ySlice[k] = (phiK %*% theta)[1]
    uSlice[k] = u[k]
    eSlice[k] = 0
  }

  return(ySlice[p:N])
}

kStepAhead = function (model, y, u, K) {
  cat('K-steap-ahead is looking into the future!!!')
}
