open_bgen <- function(bgenFile, asSampleFile = FALSE) {
  ret = .Call('_bgenR_open_bgen', PACKAGE = 'bgenR', bgenFile)
  
  if(asSampleFile) {
     if(length(ret$SampleID) == 0) {
        stop("\nERROR: There are no sample identifiers in the BGEN file.\n\n")
     } else{
        ret$SampleID <- data.frame(ID_1 = c(0, ret$SampleIdentifier),
                                   ID_2 = c(0, ret$SampleIdentifier))
     }
  }
  
  return(ret)
}



query_bgen <- function() {
  .Call('_bgenR_query_bgen', PACKAGE = 'bgenR')
}


close_bgen <- function() {
  .Call('_bgenR_close_bgen', PACKAGE = 'bgenR')
}


get_vblock <- function(){
  .Call('_bgenR_get_vblock', PACKAGE = 'bgenR')
}