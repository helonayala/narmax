
rm(list=ls()) # remove all vars and functions
cat("\014")   # clean console
while(!is.null(dev.list())) dev.off()     # clear all graphs

source("library_sysid.R")

regMatNARX = function(u,y,nu,ny,l){
  n = nu+ny
  p = max(nu,ny)+1
  auxexp = list()
  candlist = list()
  # generate all terms product combinations possible
  for(i in 1:l){
    eval(parse(text=paste0("auxexp$x",i,"=1:n"))) # generate input args for expand.grid
    cand = do.call(expand.grid,auxexp) # call expand.grid for changing number of arguments
    
    cand = t(apply(cand,1,sort)) # order each row of the matrix
    cand = unique(cand) # keep unique rows
    candlist[[i]] = cand
  }
  
  P0 = regMatrix_ARX(y,u,ny,nu,p) # creates the initial regression matrix (with ARX terms)
  P0[,1:ny] = -P0[,1:ny] # fix -y(k-..) from arx modeling
  colnames(P0) = c(paste0("y(k-",1:ny,")"),paste0("u(k-",1:nu,")")) # fix column names
  NP = nrow(P0)
  P = NULL # here we store the final NARX regression matrix 
  P = cbind(rep(1,NP), P) # constant term
  colnames(P) = "constant"
  P = cbind(P, P0) # l = 1, trivial

  for (i in 2:l){
    ncand = nrow(candlist[[i]])
    for(j in 1:ncand){
      Pcand_a = P0[,candlist[[i]][j,]]
      Pcand_b = matrix(apply(Pcand_a,1,prod),ncol = 1)
      colnames(Pcand_b) = str_c(colnames(P0[,candlist[[i]][j,]]),collapse = "")
      P = cbind(P,Pcand_b)
    }
  }
  return(P)
}

# generate simulation data ------------------------------------------------
N = 200
u = runif(N,min = -1,max = 1)
e = rnorm(N,mean = 0,sd=0.1)
y = rep(0,length(u))
for (k in 3:N) {
  y[k] = -0.605*y[k-1] - 0.163*y[k-2]^2 + 0.588*u[k-1] - 0.24*u[k-2] + e[k]
}
nu = 2
ny = 3
ne = 0
l = 3
n = nu + ny + ne

P = regMatNARX(u,y,nu,ny,l)


