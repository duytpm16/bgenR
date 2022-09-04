results <- system.file("extdata", "v1000.rds", package = "bgenR")
expected_results <- readRDS(results)

# BGENv1.1 - Uncompressed
bgen11_uncompressed <- system.file("extdata", "bgen11_uncompressed.bgen", package = "bgenR")
test_that("BGENv1.1 Uncompressed Read", {
  bgen <- open_bgen(bgen11_uncompressed)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 0)
  for (i in 1:bgen$M) { results <- query_bgen(bgen) }
  close_bgen(bgen)

  expect_equal(results, expected_results)
})


test_that("BGENv1.1 Uncompressed Seek", {
  bgen <- open_bgen(bgen11_uncompressed, getIndices = T)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 0)
  results <- query_bgen(bgen, seek = bgen$M)
  close_bgen(bgen)

  expect_equal(results, expected_results)
})


test_that("BGENv1.1 Uncompressed EOF Error", {
  bgen <- open_bgen(bgen11_uncompressed)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 0)
  expect_error(
    for (i in 1:(bgen$M + 1)) {
      results <- query_bgen(bgen)
    }
  )
  close_bgen(bgen)
})


test_that("BGENv1.1 Uncompressed EOF Error Seek", {
  bgen <- open_bgen(bgen11_uncompressed, getIndices = T)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 0)
  results <- query_bgen(bgen, seek = bgen$M)
  expect_error(query_bgen(bgen))
  close_bgen(bgen)
})



# BGENv1.1 - Zlib
bgen11_zlib <- system.file("extdata", "bgen11_zlib.bgen", package = "bgenR")
test_that("BGENv1.1 Zlib Read", {
  bgen <- open_bgen(bgen11_zlib)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 1)
  for (i in 1:bgen$M) { results <- query_bgen(bgen) }
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.1 Zlib Seek", {
  bgen    <- open_bgen(bgen11_zlib, getIndices = T)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 1)
  results <- query_bgen(bgen, seek = bgen$M)
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.1 Zlib EOF Error", {
  bgen    <- open_bgen(bgen11_zlib)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 1)
  expect_error(
    for (i in 1:(bgen$M + 1)) {
      results <- query_bgen(bgen)
    }
  )
  close_bgen(bgen)
})


test_that("BGENv1.1 Zlib EOF Error Seek", {
  bgen    <- open_bgen(bgen11_zlib, getIndices = T)
  expect_equal(bgen$layout_flag, 1)
  expect_equal(bgen$compression_flag, 1)
  results <- query_bgen(bgen, seek = bgen$M)
  expect_error(query_bgen(bgen))
  close_bgen(bgen)
})