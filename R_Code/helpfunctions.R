#Help FUnctios
library("dplyr")
library("ggplot2")
"
data        mcmc dataset (coda)
columns     specify all collumns that should get extracted of that object

returns:    specific columns of the mcmc-object in matrix form (3 lists are joined if 3 chains where initialized now)
"
get_values <- function(data,column){
  for ( i in 1:length(data)){
    data[[i]] <- data[[i]][,column]
  }
  return(as.matrix(data))
}

"
variable    which variable is planned to be subsettet (which of the collumns) (works only for one collumn or one group)
data        mcmc dataset (coda)
function_   what is suposed to happen with the subsetted collumn default is mean
...         Further options for the function_


"

get_values2 <- function(variable, data, function_ = mean,...){
  subset_pred <- grepl(paste0("^",variable,"(\\[.*\\]){0,1}$"), dimnames(data[[1]])[[2]],perl = T)
  subset_val <- get_values(data,subset_pred)
  if(length(colnames(subset_val)) == 1){colnames(subset_val) <-variable}
  return_stuff <- apply(subset_val,2,function_,...)
  return(return_stuff)
}



"
params is a list of parameters that are calculated in the chain (of the data)
where you want to extract the mean or smth similar
data: is a coda object
functions_: is a function that indicate what measure should get extracted
...         Further options for the function_

return      just makes it possible to use get_values2 multivariate
"
get_params <- function(params, data, functions_ = mean){
  a <- NULL
  for (i in 1:length(params)){
    function_results <- get_values2(params[i], data, functions_)
    #no better idea right now than this bad way o doing it...
    if(!is.null(dim(function_results))){
      a <- cbind(a,function_results)
    }else{
      a <- c(a,function_results)
    }
  }
  return((a))
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

"
makeQ creates penalization matrix for B-Splines:
https://bragqut.github.io/2016/05/24/samclifford-splines/


"
makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
}
makeQ(3,4)

