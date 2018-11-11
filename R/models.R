#' @title ARX Model
#' @description Creates an autoregressive with exogenous inputs model
#' @param ny Number of autoregressive lags
#' @param nu Number of input lags
#' @return Object representing an ARX model
#' @export
arx = function (ny, nu) {
  model = list(
    ny = ny,
    nu = nu,
    maxLag = max(ny, nu) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'arx'
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
  class(model) = 'armax'
  return(model)
}

#' @title ARMA Model
#' @description Creates an autoregressive moving average model
#' @param ny Number of autoregressive lags
#' @param ne Number of moving average lags
#' @return Object representing an ARMA model
#' @export
armax = function (ny, ne) {
  model = list(
    ny = ny,
    ne = ne,
    maxLag = max(ny, nu, ne) + 1,
    terms = NULL,
    coefficients = NULL, # Variable name required by base::
    call = match.call()
  )
  class(model) = 'arma'
  return(model)
}

#' @title NARX Model
#' @description Creates a nonlinear autoregressive with exogenous inputs model
#' @param ny Number of autoregressive lags
#' @param nu Number of input lags
#' @param nl Nonlinearity polinomial length
#' @return Object representing a NARX model
#' @export
narx = function (ny, nu, nl) {
  model = list(
    ny = ny,
    nu = nu,
    nl = nl,
    maxLag = max(ny, nu) + 1,
    terms = NULL,
    coefficients = NULL,
    call = match.call()
  )
  class(model) = 'narx'
  return(model)
}

#' @title NARMAX Model
#' @description Creates a nonlinear autogressive moving average with exogenous inputs model
#' @param ny Number of autoregressive lags
#' @param nu Number of input lags
#' @param ne Number of moving average lags
#' @param nl Nonlinearity polinomial length
#' @return Object representing a NARMAX model
#' @export
narmax = function (ny, nu, ne, nl) {
  model = list(
    ny = ny,
    nu = nu,
    ne = ne,
    nl = nl,
    maxLag = max(ny, nu, ne) + 1,
    terms = NULL,
    coefficients = NULL,
    call = match.call()
  )
  class(model) = 'narmax'
  return(model)
}

#' @title NARMA Model
#' @description Creates a nonlinear autogressive moving average model
#' @param ny Number of autoregressive lags
#' @param ne Number of moving average lags
#' @param nl Nonlinearity polinomial length
#' @return Object representing a NARMA model
#' @export
narmax = function (ny, ne, nl) {
  model = list(
    ny = ny,
    ne = ne,
    nl = nl,
    maxLag = max(ny, ne) + 1,
    terms = NULL,
    coefficients = NULL,
    call = match.call()
  )
  class(model) = 'narma'
  return(model)
}
