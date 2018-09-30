#' @title ARX Model
#' @description Creates an autoregressive with exogenous inputs model
#' @param ny Number of autoregressive lags
#' @param nu Number of input lags
#' @return Object representing an ARX model
#' @export
arx = function (ny, nu, p) {
  model = list(
    ny = ny,
    nu = nu,
    maxLag = max(ny, nu) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'arxModel'
  return(model)
}

#' @title ARMAX Model
#' @description Creates an autoregressive moving average with exogenous inputs model
#' @param ny Number of autoregressive lags
#' @param nu Number of input lags
#' @param ne Number of moving average lags
#' @return Object representing an ARMAX model
#' @export
armax = function (ny, nu, ne) {
  model = list(
    ny = ny,
    nu = nu,
    ne = ne,
    maxLag = max(ny, nu, ne) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'armaxModel'
  return(model)
}
