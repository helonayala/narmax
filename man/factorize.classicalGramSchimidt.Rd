% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/factorizers.R
\name{factorize.classicalGramSchimidt}
\alias{factorize.classicalGramSchimidt}
\title{Classical Gram-Schimidt Factorization}
\usage{
factorize.classicalGramSchimidt(X)
}
\arguments{
\item{X}{The matrix to be factorized}
}
\value{
list(Q, A) An object containing factorization results
}
\description{
Obtains \eqn{X = Q * A}, where \eqn{X in N x Nth},
so that \eqn{Q} is a \eqn{N x Nth} matrix with ortoghonal columns
and \eqn{A} is a \eqn{Nth x Nth} unit upper triangular matrix.

Reference: Aguirre 2015 Book
}
\author{
Helon Vicente Hultmann Ayala
}
