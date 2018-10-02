#' @export
predict = function (model, ...) UseMethod('predict')

#' @export
predict.arx = function (model, ...) {
  cat('Using arx prediction')
}

#' @export
predict.armax = function () {
  cat('Using armax prediction')
}

#' @export
predict.narx = function (model, ...) {
  cat('Using narx prediction')
}

#' @export
predict.narmax = function (model, y, u, K = 1, ...) {
  cat('Using narmax prediction\n')
  return(predict.default(model, y, u, K))
}

#' @export
predict.default = function (model, y, u, K, ...) {
  cat('Using default prediction\n')
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
  e = rep(0, length(y))
  P = genRegMatrix(model, y, u, e)$P
  theta = as.matrix(model$coefficients)
  print(theta)
  Yp = P %*% theta

  if (any(grepl('e(', colnames(P), fixed = TRUE))) {
    target = genTarget(model, y)
    for (i in 1:model$maxLag) {
      e = c(rep(0, model$maxLag - 1), target - Yp)
      P = genRegMatrix(model, y, u, e)$P
      Yp = P %*% theta
    }
  }

  return(Yp)
}

freeRun = function (model, y, u, ...) {

}

kStepAhead = function (model, y, u, K) {

}
