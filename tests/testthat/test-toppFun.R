###test for api entrez result
test_that("entrez lookup works", {
  expect_equal(get_Entrez("FLDB"), 338)
})

##test inputs
test_that("degs dataframe as input works",{
  data("ifnb.de")
  toppData <- toppFun(ifnb.de,
                      type="degs",
                      topp_categories="GeneOntologyMolecularFunction",
                      cluster_col="celltype",
                      gene_col="gene",
                      logFC_col = "avg_log2FC",
                      p_val_col = "p_val_adj",
                      verbose = F)
  expect_gt(nrow(toppData), 0)
})

test_that("marker df as input works",{
  data("ifnb.markers.df")
  toppData <- toppFun(ifnb.markers.df,
                      type = "marker_df",
                      topp_categories="GeneOntologyMolecularFunction",
                      verbose = F)
  expect_gt(nrow(toppData), 0)
})

test_that("marker list as input works",{
  data("ifnb.markers.list.CD8T")
  toppData <- toppFun(ifnb.markers.list.CD8T,
                      type = "marker_list",
                      topp_categories="GeneOntologyMolecularFunction",
                      verbose = F)
  expect_gt(nrow(toppData), 0)
})
##test plots
