test_that("checking estimatecc output (without .find_dmrs)", {
  library(FlowSorted.Blood.450k)
  rgset <- FlowSorted.Blood.450k[,
              pData(FlowSorted.Blood.450k)$CellTypeLong %in% "Whole blood"]
  set.seed(1234)
  est <- estimatecc(object = rgset) 
  
  expect_equal(cell_counts(est)[1,1], 0.39848, tolerance = 0.0001)
})
