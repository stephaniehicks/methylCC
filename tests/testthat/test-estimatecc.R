test_that("checking estimatecc output (without .find_dmrs)", {
  # This is a reduced version of the FlowSorted.Blood.450k 
  # dataset available by using BiocManager::install("FlowSorted.Blood.450k),
  # but for purposes of the example, we use the smaller version. 
  
  dir <- system.file("data", package="methylCC")
  files <- file.path(dir, "FlowSorted.Blood.450k.sub.RData") 
  if(file.exists(files)){
      load(file = files)
      
      set.seed(12345)
      est <- estimatecc(object = FlowSorted.Blood.450k.sub) 
      cell_counts(est)
      expect_equal(cell_counts(est)[1,1], 0.49418, tolerance = 0.0001)
  }   
 
})
