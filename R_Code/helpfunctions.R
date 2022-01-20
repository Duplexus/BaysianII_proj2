#Help FUnctios


"
data        mcmc dataset (coda)
columns       specify all collumns that should get extracted of that object
"
get_values <- function(data,column,...){
  for ( i in 1:length(data)){
    data[[i]] <- data[[i]][,column]
  }
  return(as.matrix(data)[...,])
}
"
variable    which variable is planned to be subsettet (which of the collumns)
data        mcmc dataset (coda)
function_   what is suposed to happen with the subsetted collumn default is mean

"

get_values2 <- function(variable, data, function_ = mean){
  subset_pred <- grepl(variable, dimnames(data[[1]])[[2]])
  subset_val <- get_values(data,subset_pred)
  if(length(colnames(subset_val)) == 1){colnames(subset_val) <-variable}
  return_stuff <- apply(subset_val,2,function_)
  return(return_stuff)
}

"
params is a list of parameters that are calculated in the chain (of the data)
where you want to extract the mean or smth similar
data: is a coda object
functions_: is a function that indicate what measure should get extracted
"
get_params <- function(params, data, functions_ = mean){
  a <- c()
  for (i in 1:length(params)){
    function_results <- get_values2(params[i], data, functions_)
    a <- c(a,function_results)
  }
  return(unlist(a))
}

"
variable    which variable is planned to be subsettet (which of the collumns)
data        mcmc dataset (coda)
"
subset_coda_params <- function(variable, data){
  sub_logical <- rep(F, times = length(dimnames(data[[1]])[[2]]))
  for (i in 1:length(variable)){
  subset_pred <- grepl(variable[i], dimnames(data[[1]])[[2]])
  sub_logical <- sub_logical | subset_pred
  }
  for (j in 1:length(data)){
    data[[j]]<- data[[j]][,sub_logical]
  }

  return(data)
}

"
data  output of coda
!!ATTENTION it is assumed WAIC is calculated, else you need to modify
"

summarise_default <- function(data){
  plot(data)
  cat("\n That is the effective sample size \n\n")
  print(effectiveSize(data))
  #Tine-series SE is monte carlo standard error var von dem ding durch samplsize
  cat("\n Summary of the values \n")
  print(summary(data))
  cat("\n Geweke diag and plot \n")
  print(geweke.diag(data))
  geweke.plot(data)
  cat("\n Gelman diag and plot \n")
  print(gelman.diag(data, confidence = 0.95))
  gelman.plot(data, confidence = 0.95)
  cat("\n Autocorrelation of the variables \n")
  autocorr.diag(data)
}



