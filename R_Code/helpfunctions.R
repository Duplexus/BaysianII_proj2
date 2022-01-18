#Help FUnctios

get_values <- function(object,column,...){
  for ( i in 1:length(object)){
    object[[i]] <- object[[i]][,column]
  }
  return(as.matrix(object)[...,])
}
