#' @title Clear Workspace
#'
#' @description Clear the workspace (Actually just a rdoxygen example)
#' @export
clearWorkspace = function () {
  rm(list=ls())
  cat('\014')
  while (!is.null(dev.list())) dev.off()
}

#' @title Say Hello
#'
#' @description Prints hello to the console
#' @param name Name of the user to greet
#' @example examples/sayHello.R
#' @export
sayHello = function (name, ...) {
  print(sprintf('Hello, %s!', name))
}
