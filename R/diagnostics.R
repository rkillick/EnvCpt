aic.envcpt = function(object,...,k=2){
  if(!is.list(object)){stop("object argument must be a list")}
  if(!is.matrix(object[[1]])){stop("first element in the object list must be a matrix.")}
  if(any(!is.numeric(object[[1]][c(1:2),]),na.rm=TRUE)){stop("First two rows in matrix in first element of object list must be numeric")}
  ##Calculate AIC for a secure list object
  return(object[[1]][1,] + k*object[[1]][2,])
}