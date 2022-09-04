results <- system.file("extdata", "v1000.rds", package = "bgenR")
expected_results <- readRDS(results)
expected_results$Probabilities <- expected_results$Probabilities[,-3]

# BGEN v1.2 Uncompressed
bgen12_uncompressed <- system.file("extdata", "bgen12_uncompressed.bgen", package = "bgenR")
test_that("BGENv1.2 Uncompressed Read", {
  bgen <- open_bgen(bgen12_uncompressed)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 0)
  for (i in 1:bgen$M) { results <- query_bgen(bgen) }
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Uncompressed Seek", {
  bgen <- open_bgen(bgen12_uncompressed, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 0)
  results <- query_bgen(bgen, seek = bgen$M)
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Uncompressed EOF Error", {
  bgen <- open_bgen(bgen12_uncompressed)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 0)
  expect_error(
    for (i in 1:(bgen$M + 1)) {
      results <- query_bgen(bgen)
    }
  )
  close_bgen(bgen)
})


test_that("BGENv1.2 Uncompressed EOF Error Seek", {
  bgen <- open_bgen(bgen12_uncompressed, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 0)
  results <- query_bgen(bgen, seek = bgen$M)
  expect_error(query_bgen(bgen))
  close_bgen(bgen)
})





# BGENv1.2 - Zlib
bgen12_zlib <- system.file("extdata", "bgen12_zlib.bgen", package = "bgenR")
test_that("BGENv1.2 Zlib Read", {
  bgen <- open_bgen(bgen12_zlib)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 1)
  for (i in 1:bgen$M) { results <- query_bgen(bgen) }
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Zlib Seek", {
  bgen <- open_bgen(bgen12_zlib, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 1)
  results <- query_bgen(bgen, seek = bgen$M)
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Zlib EOF Error", {
  bgen <- open_bgen(bgen12_zlib)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 1)
  expect_error(
    for (i in 1:(bgen$M + 1)) {
      results <- query_bgen(bgen)
    }
  )
  close_bgen(bgen)
})


test_that("BGENv1.2 Zlib EOF Error Seek", {
  bgen <- open_bgen(bgen12_zlib, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 1)
  results <- query_bgen(bgen, seek = bgen$M)
  expect_error(query_bgen(bgen))
  close_bgen(bgen)
})





# BGENv1.2 - Zlib
bgen12_zstd <- system.file("extdata", "bgen12_zstd.bgen", package = "bgenR")
test_that("BGENv1.2 Zstd Read", {
  bgen <- open_bgen(bgen12_zstd)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 2)
  for (i in 1:bgen$M) { results <- query_bgen(bgen) }
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Zstd Seek", {
  bgen <- open_bgen(bgen12_zstd, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 2)
  results <- query_bgen(bgen, seek = bgen$M)
  close_bgen(bgen)
  
  expect_equal(results, expected_results)
})


test_that("BGENv1.2 Zstd EOF Error", {
  bgen <- open_bgen(bgen12_zstd)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 2)
  expect_error(
    for (i in 1:(bgen$M + 1)) {
      results <- query_bgen(bgen)
    }
  )
  close_bgen(bgen)
})


test_that("BGENv1.2 Zstd EOF Error Seek", {
  bgen <- open_bgen(bgen12_zstd, getIndices = T)
  expect_equal(bgen$layout_flag, 2)
  expect_equal(bgen$compression_flag, 2)
  results <- query_bgen(bgen, seek = bgen$M)
  expect_error(query_bgen(bgen))
  close_bgen(bgen)
})