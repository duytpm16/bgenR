open_bgen <- function(bgenFile, getIndices = FALSE) {
  ret = .Call('_bgenR_open_bgen', PACKAGE = 'bgenR', bgenFile, getIndices)
  
  class(ret) <- "bgenRClass"
  return(ret)
}



query_bgen <- function(bgenList, seek = NA) {
  if(class(bgenList) != "bgenRClass") {
    stop("bgenList input must be an list object of class bgenRClass")
  }
  if(!is.list(bgenList)) {
    stop("bgenList must be a list.")    
  }
  
  
  if (is.na(seek)) {
    seek = 0
  } else {
    if ((seek < 1) || (seek > bgenList$M) || (seek %% 1 != 0)) {
      stop("seek must be a positive integer between [1, M], where M is the number of variants in the BGEN file.")
    }
  }
  
  .Call('_bgenR_query_bgen', PACKAGE = 'bgenR', bgenList, seek)
}


close_bgen <- function(bgenList) {
  if(class(bgenList) != "bgenRClass") {
    stop("Input must be an object list of class bgenRClass")
  }
  if(!is.list(bgenList)) {
    stop("bgenList must be a list.")    
  }
  
  .Call('_bgenR_close_bgen', PACKAGE = 'bgenR', bgenList)
}


variant_block <- function(bgenList){
  if(class(bgenList) != "bgenRClass") {
    stop("Input must be an object list of class bgenRClass")
  }
  if(!is.list(bgenList)) {
    stop("bgenList must be a list.")    
  }
  
  .Call('_bgenR_variant_block', PACKAGE = 'bgenR', bgenList)
}