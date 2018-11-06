#' @export
predict = function (model, ...) UseMethod('predict')

#' @export
predict.arx = function (model, y, u, K = 1, ...) {
  cat('Running arx prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.armax = function (model, y, u, K = 1, ...) {
  cat('Running armax prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.narx = function (model, y, u, K = 1, ...) {
  cat('Running narx prediction ... ')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.narmax = function (model, y, u, K = 1, ...) {
  cat('Running narmax prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.ann  = function (model, y, u, K = 1, ...) {
  cat('Running ann prediction ... ')
  prediction = predict.default.ann(model, y, u, K)
  cat('Done\n')
  return(prediction)
}

#' @export
predict.caret  = function (model, y, u, K = 1, ...) {
  cat('Running caret prediction ... ')
  prediction = predict.default.caret(model, y, u, K)
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

#' @export
predict.default.ann = function (model, y, u, K, ...) {
  if (K < 0) stop('K must be greater or equal to zero')
  method = switch(
    as.character(K),
    "1" = oneStepAhead.ann,
    "0" = freeRun.ann,
    kStepAhead.ann
  )
  return(method(model, y, u, K))
}

#' @export
predict.default.caret = function (model, y, u, K, ...) {
  if (K < 0) stop('K must be greater or equal to zero')
  method = switch(
    as.character(K),
    "1" = oneStepAhead.caret,
    "0" = freeRun.caret,
    kStepAhead.caret
  )
  return(method(model, y, u, K))
}
oneStepAhead = function (model, y, u, ...) {
  theta = as.matrix(model$coefficients)
  e = rep(0, length(y))
  p = model$maxLag
  N = length(y)

  # If e[k] does not exist on model, return the prediction
  if (!any(grepl('e(', model$terms, fixed = TRUE))) {
    P = genRegMatrix(model, y, u, e)$P
    yp = P %*% theta
    return(yp[,])
  }
  else {
    ySlice = y[1:(p-1)]
    eSlice = rep(0,p-1)

    N = length(y)

    pb = progress::progress_bar$new(total = N-p+1)
    for (k in p:N) {
      pb$tick()

      auxY = c(y     [(k - p + 1):(k - 1)], 0)
      auxU = c(u     [(k - p + 1):(k - 1)], 0)
      auxE = c(eSlice[(k - p + 1):(k - 1)], 0)
      phiK = genRegMatrix(model, auxY, auxU, auxE)$P
      ySlice[k] = (phiK %*% theta)[1]
      eSlice[k] = y[k] - ySlice[k]
    }
  }

  return(ySlice[p:N])
}

freeRun = function (model, y, u, K, ...) {
  theta = as.matrix(model$coefficients)
  p = model$maxLag
  e = rep(0, length(y))

  ySlice = y[1:(p - 1)]
  uSlice = u[1:(p - 1)]
  eSlice = e[1:(p - 1)]

  N = length(y)

  pb = progress::progress_bar$new(total = N-p+1)
  for (k in p:N) {
    pb$tick()

    auxY = c(ySlice[(k - p + 1):(k - 1)], 0)
    auxU = c(uSlice[(k - p + 1):(k - 1)], 0)
    auxE = c(eSlice[(k - p + 1):(k - 1)], 0)
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

oneStepAhead.ann = function (model, y, u, ...) {

  osa_inp = genRegMatrix(model,y,u)$P

  yh = keras::predict_on_batch(model$mdl, x = osa_inp)

  return(yh[,])
}

freeRun.ann = function (model, y, u, K, ...) {

  p = model$maxLag

  ySlice = y[1:(p - 1)]
  uSlice = u[1:(p - 1)]

  N = length(y)

  pb = progress::progress_bar$new(total = N-p+1)
  for (k in p:N) {
    pb$tick()

    auxY = c(ySlice[(k - p + 1):(k - 1)], 0)
    auxU = c(uSlice[(k - p + 1):(k - 1)], 0)
    fr_input = genRegMatrix(model, auxY, auxU)$P
    ySlice[k] = keras::predict_on_batch(model$mdl, x = t(fr_input))
    uSlice[k] = u[k]
  }

  return(ySlice[p:N])
}


kStepAhead.ann = function (model, y, u, K) {
  cat('K-steap-ahead is looking into the future!!!')
}

oneStepAhead.caret = function (model, y, u, ...) {

  osa_inp = data.frame(genRegMatrix(model,y,u)$P)

  yh = caret::predict.train(model$mdl, newdata =osa_inp)

  return(yh)
}

freeRun.caret = function (model, y, u, K, ...) {

  p = model$maxLag

  ySlice = y[1:(p - 1)]
  uSlice = u[1:(p - 1)]

  N = length(y)

  pb = progress::progress_bar$new(total = N-p+1)
  for (k in p:N) {
    pb$tick()

    auxY = c(ySlice[(k - p + 1):(k - 1)], 0)
    auxU = c(uSlice[(k - p + 1):(k - 1)], 0)
    fr_input = data.frame(t(genRegMatrix(model, auxY, auxU)$P))
    ySlice[k] = caret::predict.train(model$mdl, newdata = fr_input)
    uSlice[k] = u[k]
  }

  return(ySlice[p:N])
}

kStepAhead.caret = function (model, y, u, K) {
  cat('K-steap-ahead is looking into the future!!!')
}



