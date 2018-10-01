#' @title FROLS Algorithm
#' @description Forward Regression with Orthogonal Least Squares
#' @param P Matrix where each column is a potential feature
#' @param Y The target vector
#' @param rho Stop criteria (run while (1 - sum(EER)) > rho)
#' @export
frols = function (P, Y, rho) {
  M = ncol(P)
  NP = nrow(P)

  # 1st step ----------------------------------------------------------------
  sig = Y[, ] %*% Y[, ]
  selectTerms = NULL
  ERRvec = NULL
  gvec = NULL

  Qs = P
  g = rep(0, M)
  ERR = rep(0, M)
  for (m in 1:M) {
    g[m] = (Y[, ] %*% Qs[, m]) / (Qs[, m] %*% Qs[, m])
    ERR[m] = (g[m] ^ 2 * (Qs[, m] %*% Qs[, m])) / sig
  }
  l1 = which(ERR == max(ERR))
  selectTerms = l1[1] # vector keeping all selected terms

  # init
  # A = diag(1,M)
  A = 1
  Qs = matrix(P[, l1], ncol = 1)
  gvec = g[l1]
  ERRvec = ERR[l1]

  # s-th step --------------------------------------------------------------
  for (s in 2:M) {
    gm  = rep(0, M)
    Qm  = matrix(0, NP, M)
    ERR = rep(0, M)
    A = cbind(rbind(A, 0), 0)

    for (m in (1:M)[-selectTerms]) {
      sumQm = rep(0, NP)
      for (r in 1:(s - 1)){
        sumQm = sumQm + ((P[, m] %*% Qs[, r]) /  (Qs[, r] %*% Qs[, r])) %*% Qs[, r]
      }
      Qm[, m] = P[, m] - sumQm
      gm[m] = (Y[, ] %*% Qm[, m]) / (Qm[, m] %*% Qm[, m])
      ERR[m] = (gm[m] ^ 2 * (Qm[, m] %*% Qm[, m])) / sig
    }

    ls = which(ERR == max(ERR))
    selectTerms = cbind(selectTerms, ls) # vector keeping all selected terms

    Qs = cbind(Qs, Qm[, ls]) # keep set of orthogonal bases
    gvec = rbind(gvec, gm[ls])
    for (r in 1:(s - 1)){
      A[r, s] = (Qs[, r] %*% P[, ls]) / (Qs[, r] %*% Qs[, r])
    }
    A[s, s] = 1

    ERRvec = rbind(ERRvec, ERR[ls])
    ESR = 1 - sum(ERRvec)

    if (ESR <= rho){
      M0 = s
      break
    }
  }

  th_FROLS = solve(A, gvec)
  Psel = P[, selectTerms]
  return(list(
    th = th_FROLS[, ],
    Psel = Psel,
    g = gvec,
    W = Qs,
    A = A,
    ERR = ERRvec[, ]
  ))
}
