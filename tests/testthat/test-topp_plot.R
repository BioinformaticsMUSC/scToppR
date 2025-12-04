###test for toppPlot
test_that("toppPlot works",{
  data("ifnb.de")
  toppData <- toppFun(ifnb.de,
                      type="degs",
                      topp_categories="GeneOntologyMolecularFunction",
                      cluster_col="celltype",
                      gene_col="gene",
                      logFC_col = "avg_log2FC",
                      p_val_col = "p_val_adj",
                      verbose = FALSE)

    dotplot <-toppPlot(toppData,
             category="GeneOntologyMolecularFunction",
             clusters=0,
             save=FALSE)

    expect_s3_class(dotplot, "ggplot")
    expect_equal(length(dotplot$layers), 2)

})

test_that("toppPlot multiple clusters works",{
  data("ifnb.de")
  toppData <- toppFun(ifnb.de,
                      type="degs",
                      topp_categories="GeneOntologyMolecularFunction",
                      cluster_col="celltype",
                      gene_col="gene",
                      logFC_col = "avg_log2FC",
                      p_val_col = "p_val_adj",
                      verbose = FALSE)

    dotplot_list <-toppPlot(toppData,
             category="GeneOntologyMolecularFunction",
             clusters=c("CD4T","CD8T"),
             save=FALSE)

    expect_type(dotplot_list, "list")
    expect_equal(length(dotplot_list), 2)
    expect_s3_class(dotplot_list$CD4T, "ggplot")
    expect_equal(length(dotplot_list$CD4T$layers), 2)
})

test_that("toppBalloon works",{
  data("ifnb.de")
  toppData <- toppFun(ifnb.de,
                      type="degs",
                      topp_categories="GeneOntologyMolecularFunction",
                      cluster_col="celltype",
                      gene_col="gene",
                      logFC_col = "avg_log2FC",
                      p_val_col = "p_val_adj",
                      verbose = FALSE)

    balloonplot <-toppBalloon(toppData,
             categories="GeneOntologyMolecularFunction",
             save=FALSE)

    expect_s3_class(balloonplot, "ggplot")
    expect_equal(length(balloonplot$layers), 1)

})