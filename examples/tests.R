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

A = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
colnames(A) = c('ColA', 'ColB', 'ColC')
rownames(A) = c('k=1', 'k=2', 'k=3')

print(A)
