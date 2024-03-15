test_that("RcppKappaCluster returns expected output", {
  # Read the data files
  gs1 <- read.delim(system.file("extdata", "go1.txt", package = "RichStudio"))
  gs2 <- read.delim(system.file("extdata", "go2.txt", package = "RichStudio"))
  gs3 <- read.delim(system.file("extdata", "kegg1.txt", package = "RichStudio"))

  # Merge the gene sets
  mergelist <- list(gs1, gs2, gs3)
  names(mergelist) <- c("gs1", "gs2", "kegg1")
  merge <- RichStudio::merge_genesets(mergelist)

  # Call the RcppKappaCluster function
  result <- RichStudio::RcppKappaCluster(merge$Term, merge$GeneID)

  # Test the second version of it!
  result2 <- RichStudio::RcppKappaCluster2(merge$Term, merge$GeneID, merge$Padj)

  # View files in window
  View(result[["SignifKappaTerms"]])
  View(result[["SignifKappaScores"]])
  View(result[["KappaScores"]])

  View(result2[["dfSigTermIndices"]])
  View(result2[["all_kappas"]])
  View(result2[["AllKappaMatrix"]])

  # Compare to richR
  richRresult <- richR::richCluster(x=merge, gene=TRUE, cutoff=0.5, overlap=0.5, minSize=2)

  richr_sigterms <- merge[merge$Annot %in% c("GO:0043624","GO:1902904","GO:0051261","GO:0032984","GO:0043241"),]

  # Check the structure of the result
  expect_is(result, "list")
  expect_named(result, c("SignifKappaTerms", "KappaScores"))
  expect_s3_class(result$SignifKappaTerms, "data.frame")
  expect_s3_class(result$KappaScores, "data.frame")

  # You can add more specific checks for the content of the data frames if needed
})
