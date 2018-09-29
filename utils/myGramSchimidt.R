# Classical Gram-Schmidt factorization
# helon - 4/9/18
# mec2015 - system identification - puc-rio

# obtains Phi = Q * A
# Phi is N x Nth
# so that
# Q is a N x Nth matrix with orthogonal columns 
# A is a Nth x Nth unit upper triangular matrix

Mat = matrix(c(9,-5,5,-1.53,9.08,
               1,-1,8,-0.39,7.87,
               2,-5,6,-3.26,3.01,
               8,-2,0,0.36,5.98,
               0,0,9,0.13,9.05), # the data elements 
             nrow=5,              # number of rows 
             byrow = TRUE)        # fill matrix by rows

P = Mat[,1:4]

N   = nrow(P)
Nth = ncol(P)

# init matrix
A = diag(1,nrow = Nth)
Q = matrix(0,nrow =N, ncol = Nth)
Q[,1] = P[,1]
for (i in 2:Nth) {
  Q[,i] = P[,i] 
  for (j in 1:(i-1)){
    print(j,i)
    A[j,i] = (Q[,j] %*% P[,i]) / (Q[,j] %*% Q[,j])
    Q[,i] = Q[,i] - A[j,i] * Q[,j]
  }
}




