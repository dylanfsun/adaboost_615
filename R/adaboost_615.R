

adaboost_615 <- function(data, dependent_index) {
  if(missing(dependent_index)){
    dependent_index <- dim(data)[2]
  } else {
    data <- data[c(1:dependent_index-1, dependent_index+1:dim(data)[2], dependent_index)]
  }


  data <- as.matrix(data)
  return(adaboost(data))
}
