test_that("RcppKappaCluster returns expected output", {
  # Read the data files
  gs1 <- read.delim(system.file("extdata", "go1.txt", package = "RichStudio"))
  gs2 <- read.delim(system.file("extdata", "go2.txt", package = "RichStudio"))

  # Merge the gene sets
  mergelist <- list(gs1, gs2)
  names(mergelist) <- c("gs1", "gs2")
  merge <- RichStudio::merge_genesets(mergelist)

  # Call the RcppKappaCluster function
  result <- RichStudio::RcppKappaCluster(merge$Term, merge$GeneID)
  richRresult <- richR::richCluster(x=merge, gene=TRUE)

  # Check the structure of the result
  expect_is(result, "list")
  expect_named(result, c("SignifKappaTerms", "KappaScores"))
  expect_s3_class(result$SignifKappaTerms, "data.frame")
  expect_s3_class(result$KappaScores, "data.frame")

  # You can add more specific checks for the content of the data frames if needed
})
