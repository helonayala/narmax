#' @export
predict = function (model, ...) UseMethod('predict')

#' @export
predict.arx = function (model, y, u, K = 1, ...) {
  cat('Running arx prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
  return(prediction)
}

#' @export
predict.armax = function (model, y, u, K = 1, ...) {
  cat('Running armax prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
  return(prediction)
}

#' @export
predict.nar = function (model, y, u, K = 1, ...) {
  cat('Running nar prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
  return(prediction)
}

#' @export
predict.narma = function (model, y, K = 1, ...) {
  cat('Running narma prediction ... \n')
  prediction = predict.default.narma(model, y, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
  return(prediction)
}

#' @export
predict.narx = function (model, y, u, K = 1, ...) {
  cat('Running narx prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
  return(prediction)
}

#' @export
predict.narmax = function (model, y, u, K = 1, ...) {
  cat('Running narmax prediction ... \n')
  prediction = predict.default(model, y, u, K)
  cat(sprintf('Done. R2 = %0.4f\n',prediction$R2))
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
predict.default.narma = function (model, y, K, ...) {
  if (K < 0) stop('K must be greater or equal to zero')
  method = switch(
    as.character(K),
    "1" = oneStepAhead.narma,
    "0" = freeRun.narma,
    kStepAhead.narma
  )
  return(method(model, y, K))
}

oneStepAhead = function (model, y, u,...) {
  theta = as.matrix(model$coefficients)
  e = rep(0, length(y))
  p = model$maxLag
  N = length(y)

  type = "one-step-ahead"
  if (class(mdl) %in% c("armax","arx")) nl = 1 # linear with X
  else if (class(mdl) %in% "narmax")    nl = 2 # nonlinear with x
  else                                  nl = 3 # undefined

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
  df = data.frame(time = p:N,
                  y = y[p:N],
                  u = u[p:N],
                  yh = ySlice[p:N],
                  e = y[p:N] - ySlice[p:N])
  g = xcorrel(df$e,df$u,nl)

  R2 = calcR2(y[p:N],ySlice[p:N])

  p = cookPlots(df,R2,type)

  out = list(dfpred = df,
             R2 = R2,
             ploty = p[[1]],
             plote = p[[2]],
             xcorrel = g,
             type = type)

  return(out)
}

freeRun = function (model, y, u, K, ...) {
  theta = as.matrix(model$coefficients)
  p = model$maxLag
  e = rep(0, length(y))
  type = "free-run"

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

  df = data.frame(time = p:N,
                  y = y[p:N],
                  u = u[p:N],
                  yh = ySlice[p:N],
                  e = y[p:N] - ySlice[p:N])

  R2 = calcR2(y[p:N],ySlice[p:N])

  p = cookPlots(df,R2,type)

  out = list(dfpred = df,
             R2 = R2,
             ploty = p[[1]],
             plote = p[[2]],
             type = type)
  return(out)
}

kStepAhead = function (model, y, u, K) {
  cat('K-steap-ahead is looking into the future!!!')
}

oneStepAhead.narma = function (model, y, ...) {
  theta = as.matrix(model$coefficients)
  e = rep(0, length(y))
  p = model$maxLag
  N = length(y)

  type = "one-step-ahead"

  # If e[k] does not exist on model, return the prediction
  if (!any(grepl('e(', model$terms, fixed = TRUE))) {
    P = genRegMatrix(model, Y = y, E = e)$P
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
      auxE = c(eSlice[(k - p + 1):(k - 1)], 0)
      phiK = genRegMatrix(model, Y = auxY, E = auxE)$P
      ySlice[k] = (phiK %*% theta)[1]
      eSlice[k] = y[k] - ySlice[k]
    }
  }

  df = data.frame(time = p:N,
                  y = y[p:N],
                  yh = ySlice[p:N],
                  e = y[p:N] - ySlice[p:N])

  g = xcorrel.ts(df$e)

  R2 = calcR2(y[p:N],ySlice[p:N])

  p = cookPlots(df,R2,type)

  out = list(dfpred = df,
             R2 = R2,
             ploty = p[[1]],
             plote = p[[2]],
             xcorrel = g,
             type = type)

  return(out)
}

freeRun.narma = function (model, y, K, ...) {
  theta = as.matrix(model$coefficients)
  p = model$maxLag
  type = "free-run"
  e = rep(0, length(y))
  ySlice = y[1:(p - 1)]
  eSlice = e[1:(p - 1)]
  N = length(y)

  pb = progress::progress_bar$new(total = N-p+1)
  for (k in p:N) {
    pb$tick()
    auxY = c(ySlice[(k - p + 1):(k - 1)], 0)
    auxE = c(eSlice[(k - p + 1):(k - 1)], 0)
    phiK = genRegMatrix(model, Y = auxY,E=auxE)$P
    ySlice[k] = (phiK %*% theta)[1]
  }

  df = data.frame(time = p:N,
                  y = y[p:N],
                  yh = ySlice[p:N],
                  e = y[p:N] - ySlice[p:N])

  R2 = calcR2(y[p:N],ySlice[p:N])

  p = cookPlots(df,R2,type)

  out = list(dfpred = df,
             R2 = R2,
             ploty = p[[1]],
             plote = p[[2]],
             type = type)
  return(out)
}

kStepAhead.narma = function (model, y, K) {

  theta = as.matrix(model$coefficients)
  p = model$maxLag
  N = length(y)
  e = rep(0, N)
  type = paste0(K,"-steps ahead")

  pb = progress::progress_bar$new(total = N-p-K+2)

  time = (p+K-1):N
  yh = 0*time

  for (k in p:(N-K+1)) {
    pb$tick()

    ySlice = y[(k - p + 1):(k - 1)]
    eSlice = e[(k - p + 1):(k - 1)]
    for (nsteps in 1:K){
      auxY = c(ySlice[(nsteps):(p+nsteps-2)], 0)
      auxE = c(eSlice[(nsteps):(p+nsteps-2)], 0)
      phiK = genRegMatrix(model,Y = auxY, E = auxE)$P
      ySlice[p + nsteps - 1] = (phiK %*% theta)[1]
      eSlice[p + nsteps - 1] = 0
    }
    yh[k-p+1] = ySlice[p + K - 1]
  }

  df = data.frame(time = (p+K-1):N,
                  y = y[(p+K-1):N],
                  yh = yh,
                  e = y[(p+K-1):N] - yh)

  R2 = calcR2(df$y,df$yh)

  p = cookPlots(df,R2,type)

  out = list(dfpred = df,
             R2 = R2,
             ploty = p[[1]],
             plote = p[[2]],
             type = type)

  return(out)
}
