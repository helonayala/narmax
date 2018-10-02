#' @title Clear Workspace
#'
#' @description Clear the workspace (Actually just a rdoxygen example)
#' @export
clearWorkspace = function () {
  rm(list=ls())
  cat('\014')
  while (!is.null(dev.list())) dev.off()
}

#' @title Subset matrix
#' @description Subset a matrix like mat[rows, cols]. Return a matrix even when the
#' subset generates a lower dimensional structure, i.e., never loose rownames or colnames.
#' @param mat Target matrix
#' @param rows Rows subsetter
#' @param cols Cols subsetter
#' @export
subsetMatrix = function (mat, rows, cols) {
  if (is.null(rows)) rows = 1:nrow(mat)
  if (is.null(cols)) cols = 1:ncol(mat)
  names = dimnames(mat)
  submat = matrix(
    mat[rows, cols],
    nrow = length(rows),
    ncol = length(cols),
    dimnames = list(names[[1]][rows], names[[2]][cols])
  )
  return(submat)
}

#' @title Find term indexes
#' @description Find term indexes in regression matrix
#' @param P Regression matrix
#' @param terms Vector of terms (string format)
findTermIndexes = function (P, terms) {
  names = colnames(P)
  indexes = which(names == terms[1])
  for (i in 2:length(terms)) {
    indexes = c(indexes, which(names == terms[i]))
  }
  return(indexes)
}
