allPossibleCombinations <- function(elements,repetitions=T){
#' This function generates a vector with all possible combiantions of elements
#' of a given vector of elements. Unique values will not be checked and combi-
#' nations will respect the order of elements vector. Direction doesn't matter,
#' it means that (AB it's equal tha BA)
#'  @param elements vector with elements to perform combination
#'  @param repetitions use or not use repetitions of elements (Ej: AA)
#'  @return a vector with all possible combinations or NULL if any error occurs
#' @author Fernando Moreno Jabato <jabato@uma.com>
#' @description function to obtain all possible combinations between elements of a given vector
#' @RVersion 3.4.3
  
  # Check input
  if(!is.vector(elements)){
    message("Error: elements variable must be a vector")
    return(NULL)
  }else if(!is.logical(repetitions)){
    message("Error: repetitions indicator must be logical value (TRUE,FALSE)")
    return(NULL)
  }
  
  # Obtain all possible combinations
  combinations <- expand.grid(elements,elements)
  # Alternate columns
  aux <- combinations[,2]
  combinations[,2] <- combinations[,1]
  combinations[,1] <- aux
  # Clean repetitions
  if(repetitions)
    combinations <- combinations[as.numeric(combinations[,1]) <= as.numeric(combinations[,2]),]
  else
    combinations <- combinations[as.numeric(combinations[,1]) < as.numeric(combinations[,2]),]
  
  # Reset rownames
  rownames(combinations) <- NULL
  
  # Return
  return(combinations)
}
