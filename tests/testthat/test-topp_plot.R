### test for toppPlot
test_that("toppPlot works", {
    data("toppdata.pbmc")

    dotplot <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = 0,
        save = FALSE
    )

    expect_s3_class(dotplot, "ggplot")
    expect_equal(length(dotplot$layers), 2)
})

test_that("toppPlot multiple clusters works", {
    data("toppdata.pbmc")

    dotplot_list <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = c("CD4T", "CD8T"),
        save = FALSE
    )

    expect_type(dotplot_list, "list")
    expect_equal(length(dotplot_list), 2)
    expect_s3_class(dotplot_list$CD4T, "ggplot")
    expect_equal(length(dotplot_list$CD4T$layers), 2)
})

test_that("toppBalloon works", {
    data("toppdata.pbmc")

    balloonplot <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        save = FALSE
    )

    expect_s3_class(balloonplot, "ggplot")
    expect_equal(length(balloonplot$layers), 1)
})
