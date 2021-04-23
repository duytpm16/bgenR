open_bgen <- function(bgenFile, sampleFile = FALSE, bytes = FALSE) {
  ret = .Call('_bgenR_open_bgen', PACKAGE = 'bgenR', bgenFile, bytes)
  
  if(sampleFile) {
     if(length(ret$SampleID) == 0) {
        warning("\nWARNING: There are no sample identifiers in the BGEN file.\n\n", immediate. = TRUE)
     } else{
        ret$SampleID <- data.frame(ID_1 = c(0, ret$SampleIdentifier),
                                   ID_2 = c(0, ret$SampleIdentifier))
     }
  }
  class(ret) <- "bgenRClass"
  return(ret)
}



query_bgen <- function(bgenList, seek = 0) {
  if(class(bgenList) != "bgenRClass") {
    stop("\nERROR: bgenList input must be an list object of class bgenRClass")
  }
  if(!is.list(bgenList)) {
    stop("\nERROR: bgenList must be a list.")    
  }
  
  if (seek != 0) {
    if (seek < 1 || (seek > bgenList$M)) {
      stop("\nERROR: seek must bean positive integer starting from 1 to M, where M is the number of variants in the BGEN file.")
    }
    if (seek %% 1 != 0){
      stop("\nERROR: seek input must be an integer")   
    }
  }

  
  .Call('_bgenR_query_bgen', PACKAGE = 'bgenR', bgenList, seek)
}


close_bgen <- function(bgenList) {
  if(class(bgenList) != "bgenRClass") {
    stop("\nERROR: Input must be an object list of class bgenRClass")
  }
  .Call('_bgenR_close_bgen', PACKAGE = 'bgenR', bgenList)
}


get_vblock <- function(bgenList){
  if(class(bgenList) != "bgenRClass") {
    stop("\nERROR: Input must be an object list of class bgenRClass")
  }
  
  .Call('_bgenR_get_vblock', PACKAGE = 'bgenR', bgenList)
}