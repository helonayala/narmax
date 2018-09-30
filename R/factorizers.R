#' @title Factorize matrix
#' @description Decompose a matrix \eqn{X} into two simpler matrices
#' \eqn{Q} (Orthogonal columns) and \eqn{A} (Unit upper triangular), satisfying \eqn{X = Q * A}
#' @param X The matrix to be factorized
#' @param method Which method to use. Available methods are:
#' \describe{
#'   \item{"cgs"}{Classical Gram-Schimidt}
#'   \item{"mgs"}{Modified Gram-Schimidt}
#' }
#' It uses method = "mgs" by default.
#' @return list(Q, A) Object containing factorization results
#' @export
factorize = function (X, method = "mgs") {
  return(switch (method,
    "cgs" = factorize.classicalGramSchimidt,
    "mgs" = factorize.modifiedGramSchimidt,
    factorize.modifiedGramSchimidt
  )(X))
}

#' @title Classical Gram-Schimidt Factorization
#' @description Obtains \eqn{X = Q * A}, where \eqn{X in N x Nth},
#' so that \eqn{Q} is a \eqn{N x Nth} matrix with ortoghonal columns
#' and \eqn{A} is a \eqn{Nth x Nth} unit upper triangular matrix.
#'
#' Reference: Aguirre 2015 Book
#' @author Helon Vicente Hultmann Ayala
#' @param X The matrix to be factorized
#' @export
#' @return list(Q, A) An object containing factorization results
factorize.classicalGramSchimidt = function (X) {
  P = as.matrix(X)
  N = nrow(P)
  Nth = ncol(P)

  # Initialize matrices
  A = diag(1, nrow = Nth)
  Q = matrix(0, nrow = N, ncol = Nth)
  Q[, 1] = P[, 1]

  # Factorization loop
  for (i in 2:Nth) {
    Q[, i] = P[, i]
    for (j in 1:(i - 1)) {
      A[j, i] = (Q[, j] %*% P[, i]) / (Q[, j] %*% Q[, j])
      Q[, i] = Q[, i] - A[j, i] * Q[, j]
    }
  }

  colnames(Q) = colnames(X)
  return(list(
    Q = Q,
    A = A
  ))
}

#' @title Modified Gram-Schimidt
#' @description Obtains \eqn{X = Q * A}, where \eqn{X in N x Nth},
#' so that \eqn{Q} is a \eqn{N x Nth} matrix with ortoghonal columns
#' and \eqn{A} is a \eqn{Nth x Nth} unit upper triangular matrix.
#'
#' Reference: Aguirre 2015 Book
#' @author Helon Vicente Hultmann Ayala
#' @param X The matrix to be factorized
#' @export
#' @return list(Q, A) An object containing factorization results
factorize.modifiedGramSchimidt = function (X) {
  P = as.matrix(X)
  N = nrow(P)
  Nth = ncol(P)

  # Initialize matrices
  A = diag(1, nrow = Nth)
  P_i_1 = P
  Q = matrix(0, nrow =N, ncol = Nth)
  P_i  = matrix(0, nrow = N, ncol = Nth)

  # Factorization loop
  for (i in 1:(Nth - 1)) {
    Q[,i] = P_i_1[, i]
    for (j in (i + 1):Nth){
      #disp(j,i)
      A[i, j] = (Q[, i] %*% P_i_1[, j]) / (Q[, i] %*% Q[, i])
      P_i[, j] = P_i_1[, j] - A[i, j] * Q[, i]
    }
    P_i_1 = P_i
  }

  # The end
  Q[, Nth] = P_i_1[, Nth]

  colnames(Q) = colnames(X)
  return(list(
    Q = Q,
    A = A
  ))
}
