#' @title Base Representation
#' @description Generic function to print header
#' @param model A sysident model
#' @return Boolean indicating if model is uninitialized
baseRepr = function (model) {
  print(model$call)
  hasBeenEstimated = !is.null(model$terms)
  if (!hasBeenEstimated) invisible(print("Not estimated"))
  return(hasBeenEstimated)
}

#' @title ARX model printer
#' @description Print basic info from ARX model
#' @param model ARX model
#' @export
print.arx = function (model) {
  if (baseRepr(model)) {
    invisible(print('TODO: Show terms and coefficients'))
  }
}

#' @title ARMAX model printer
#' @description Print basic info from ARMAX model
#' @param model ARMAX model
#' @export
print.armax = function (model) {
  if (baseRepr(model)) {
    invisible(print('TODO: Show terms and coefficients'))
  }
}

#' @title NARX model printer
#' @description Print basic info from NARX model
#' @param model NARX model
#' @export
print.narx = function (model) {
  if (baseRepr(model)) {
    invisible(print('TODO: Show terms and coefficients'))
  }
}

#' @title NARMAX model printer
#' @description Print basic info from NARMAX model
#' @param model NARMAX model
#' @export
print.narmax = function (model) {
  if (baseRepr(model)) {
    invisible(print('TODO: Show terms and coefficients'))
  }
}
