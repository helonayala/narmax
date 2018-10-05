
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

calc_R2 = function(y_true,y_pred){
  # calcula o R2 - coeficiente de correlação múltipla
  SSE = sum((y_true - y_pred)^2)
  avg_real = mean(y_true)
  sum2 = sum( (y_true - avg_real)^2)
  R2 = 1 - SSE / sum2;
  return(R2)
}

build_model <- function(lr) {
  model <- keras_model_sequential() %>% 
    layer_dense(units = 8, activation = "tanh", 
                input_shape = dim(train_data)[[2]]) %>% 
    layer_dense(units = 8, activation = "tanh") %>%
    layer_dense(units = 8, activation = "tanh") %>% 
    # layer_dense(units = 8, activation = "tanh") %>% 
    # layer_dense(units = 8, activation = "tanh") %>% 
    # layer_dense(units = 20, activation = "tanh") %>% 
    # layer_dense(units = 20, activation = "tanh") %>% 
    # layer_dense(units = 20, activation = "tanh") %>% 
    layer_dense(units = 1) 
  # layer_dense(units = 8, activation = "tanh", 
  #             input_shape = dim(train_data)[[2]]) %>% 
  #   layer_dense(units = 8, activation = "tanh") %>% 
  #   layer_dense(units = 8, activation = "tanh") %>% 
  #   layer_dense(units = 8, activation = "tanh") %>% 
  #   layer_dense(units = 1) 
  
  model %>% compile(
    #optimizer = optimizer_rmsprop(lr=lr), 
    optimizer = optimizer_adam(lr=lr,amsgrad = TRUE,clipnorm=50), 
    loss = "mse", 
    metrics = c("mae")
  )
}

create_reg_matrix = function(y,u,ny,nu){
  # create a matrix with specific set of regressors for each of the input/output
  # y: dados de saída
  # u: dados de entrada
  # ny: conjunto de lags em y
  # nu: conjunto de lags em u
  # X: conjunto de entradas
  # Y: saidas correspondentes a X
  
  N_tot = length(y)
  max_y_lag = max(ny)
  max_u_lag = max(nu)
  max_lag = max(c(max_y_lag,max_u_lag))
  Y = y[(max_lag+1):N_tot]
  N = length(Y)
  
  X = c()
  col_str = c()
  
  for (i in 1:max_y_lag) {
    X = cbind(X, y[(N_tot-N-(i-1)) : ( (N_tot-1)-(i-1))])
    col_str = c(col_str,paste0("y(t-",i,")"))
  }
  for (j in 1:max_u_lag) {
    X = cbind(X, u[ (N_tot-N-(j-1) ): ((N_tot-1)-(j-1))])
    col_str = c(col_str,paste0("u(t-",j,")"))
  }
  
  # atribui o nome das colunas (p/ facilitar debug)
  colnames(X) = col_str
  
  # elimina os lags não utilizados
  ind_vec = c(is.element(1:max_y_lag,ny), is.element(1:max_u_lag,nu))
  X = X[,ind_vec]
  
  return(list (X, Y) )
}

predictFreeRun = function (u,y,ny,nu,model){
  
  ndata = length(y)
  
  maxn = max(ny,nu)
  
  y2 = y[1:(maxn)]
  
  yh_fr = c()
  for (i in (maxn) : (ndata-1)) {
    list[fr_input, dummy] = 
      create_reg_matrix(c(y2,-1),u[1:i],ny,nu) #-1 é dummy
    
    fr_input2 = matrix(
      matrix(fr_input,ncol = length(ny)+length(nu))[i-maxn+1,], 
      ncol = length(ny)+length(nu)
      ) # transf. p/ matriz (1a iter) e pega a ultima linha
    
    yh_fr [i-maxn+1] = y2[i+1] = predict_on_batch(model, x = fr_input2)
    
    #print(paste("Free-run! i =",i,"de",ndata-1))
  }
  return(yh_fr)
}

