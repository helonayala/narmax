#' @title ARX Model
#' @description Creates an autoregressive with exogenous inputs model
#' @param na Number of autoregressive lags
#' @param nb Number of input lags
#' @param p Autoregressive delay ??
#' @return Object representing an ARX model
#' @export
model.arx = function (na, nb, p) {
  model = list(
    na = na,
    nb = nb,
    p = p,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'arxModel'
  return(model)
}

#' @title ARMAX Model
#' @description Creates an autoregressive moving average with exogenous inputs model
#' @param na Number of autoregressive lags
#' @param nb Number of input lags
#' @param nc Number of moving average lags
#' @param p Autoregressive delay ??
#' @return Object representing an ARMAX model
#' @export
model.armax = function (na, nb, nc, p) {
  model = list(
    na = na,
    nb = nb,
    nc = nc,
    p = p,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'armaxModel'
  return(model)
}
