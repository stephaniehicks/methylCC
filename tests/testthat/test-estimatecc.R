test_that("checking estimatecc output (without .find_dmrs)", {
  suppressPackageStartupMessages(library(FlowSorted.Blood.450k))
  data(FlowSorted.Blood.450k)
  # take a random sample to make object size in build smaller
  set.seed(12345)
  cpg_ids <- sample(seq_len(nrow(FlowSorted.Blood.450k)), 2e5)
  rgset <- FlowSorted.Blood.450k[cpg_ids,
                                 pData(FlowSorted.Blood.450k)$CellTypeLong %in% "Whole blood"]
  set.seed(12345)
  est <- estimatecc(object = rgset) 
  expect_equal(cell_counts(est)[1,1], 0.46519, tolerance = 0.0001)
})
