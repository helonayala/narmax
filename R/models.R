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

#' @title ANN Model
#' @description Creates a NARX artificial neural network nonlinear
#' model (with autogressive and exogenous inputs)
#' @param oy Vector with the lags in the autoregressive terms
#' @param ou Vector with the lags in the exogenous terms
#' @param nrn a vector (with dimension = number of hidden layers) with the number of neurons in each hidden layers
#' @param afc activation functions of the neurons (see ? layer_dense)
#' @return Object representing the ANN-NARX model
#' @export
ann = function (oy, ou, nrn, afc) {

  hdnl = length(nrn) # number of hidden layers
  ninp = length(oy)+length(ou) # number of model inputs

  mdl = keras::layer_dense(keras::keras_model_sequential() ,units = nrn[1], activation = afc,
                      input_shape = ninp, name = "hidden_layer_1")

  if (hdnl > 1){
    for (i in 2:hdnl){
      mdl =  keras::layer_dense(mdl, units = nrn[i], activation = afc, name = paste0("hidden_layer_",i))
    }
  }

  mdl = keras::layer_dense(mdl,units = 1, name = "output_layer") # 1-output layer

  model = list(
    oy = oy,
    ou = ou,
    maxLag = max(oy, ou) + 1,
    nrn = nrn,
    afc = afc,
    mdl = mdl,
    call = match.call()
  )
  class(model) = 'ann'
  return(model)
}


#' @title ANN-ts Model
#' @description Creates a ANN-ts nonlinear model (with autogressive inputs)
#' @param oy Vector with the lags in the autoregressive terms
#' @param nrn a vector (with dimension = number of hidden layers) with the number of neurons in each hidden layers
#' @param afc activation functions of the neurons (see ? layer_dense)
#' @return Object representing the ANN-ts model
#' @export
annts = function (oy, nrn, afc) {

  hdnl = length(nrn) # number of hidden layers
  ninp = length(oy)  # number of model inputs

  mdl = keras::layer_dense(keras::keras_model_sequential() ,units = nrn[1], activation = afc,
                           input_shape = ninp, name = "hidden_layer_1")

  if (hdnl > 1){
    for (i in 2:hdnl){
      mdl =  keras::layer_dense(mdl, units = nrn[i], activation = afc, name = paste0("hidden_layer_",i))
    }
  }

  mdl = keras::layer_dense(mdl,units = 1, name = "output_layer") # 1-output layer

  model = list(
    oy = oy,
    maxLag = max(oy) + 1,
    nrn = nrn,
    afc = afc,
    mdl = mdl,
    call = match.call()
  )
  class(model) = 'annts'
  return(model)
}
