read_bgen <- function(x) {
  .Call('_bgenR_set_bgen', PACKAGE = 'bgenR', x)
}


query_bgen <- function() {
  .Call('_bgenR_query_bgen', PACKAGE = 'bgenR')
}


close_bgen <- function() {
  .Call('_bgenR_close_bgen', PACKAGE = 'bgenR')
}