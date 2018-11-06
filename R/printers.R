#' @title Base Representation
#' @description Generic function to print header
#' @param model A sysident model
#' @return Boolean indicating if model is uninitialized
baseRepr = function (model) {
  print(model$call)

  ny = model$ny
  nu = model$nu
  ne = model$ne
  nl = model$nl

  props = list()
  if (!is.null(ny)) props$ny = ny
  if (!is.null(nu)) props$ny = nu
  if (!is.null(ne)) props$ny = ne
  if (!is.null(nl)) props$ny = nl

  if (is.null(model$terms)) {
    cat('\nNot estimated\n')
  } else {
    nTerms = length(model$terms)
    cat('\n')
    cat(sprintf('%20s%20s\n', 'Term', 'Coefficient'))
    for (i in 1:nTerms) {
      cat(sprintf('%20s%20.4f\n', model$terms[i], model$coefficients[i]))
    }
  }
}

#' @title ARX model printer
#' @description Print basic info from ARX model
#' @param model ARX model
#' @export
print.arx = function (model) baseRepr(model)

#' @title ARMAX model printer
#' @description Print basic info from ARMAX model
#' @param model ARMAX model
#' @export
print.armax = function (model) baseRepr(model)

#' @title NARX model printer
#' @description Print basic info from NARX model
#' @param model NARX model
#' @export
print.narx = function (model) baseRepr(model)

#' @title NARMAX model printer
#' @description Print basic info from NARMAX model
#' @param model NARMAX model
#' @export
print.narmax = function (model) baseRepr(model)

#' @title ANN-NARX model printer
#' @description Print basic info from ANN-NARX model
#' @param model ann model
#' @export
print.ann = function(model){
  print(model$mdl)
}

#' @title caret-NARX model printer
#' @description Print basic info from caret-NARX model
#' @param model caret model
#' @export
print.caret = function(model){
  print(model$mdl)
}





