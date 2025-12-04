### test for api entrez result
test_that("entrez lookup works", {
    expect_equal(get_Entrez("FLDB"), 338)
})

## test inputs
test_that("degs dataframe as input works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )
    expect_gt(nrow(toppData), 0)
})

test_that("marker df as input works", {
    data("ifnb.markers.df")
    toppData <- toppFun(ifnb.markers.df,
        type = "marker_df",
        topp_categories = "GeneOntologyMolecularFunction",
        verbose = FALSE
    )
    expect_gt(nrow(toppData), 0)
})

test_that("marker list as input works", {
    data("ifnb.markers.list.CD8T")
    toppData <- toppFun(ifnb.markers.list.CD8T,
        type = "marker_list",
        topp_categories = "GeneOntologyMolecularFunction",
        verbose = FALSE
    )
    expect_gt(nrow(toppData), 0)
})
## test save file as xlsx - split by cluster
test_that("toppFun save as split xlsx works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = TRUE,
        format = "xlsx",
        verbose = FALSE
    )
    expect_true(length(list.files(tmp_dir, pattern = "\\.xlsx$")) > 0)
    df <- openxlsx::read.xlsx(file.path(tmp_dir, "test_toppFun_CD8_T.xlsx"), sheet = "toppData")
    expect_gt(nrow(df), 0)
})

## test save file as xlsx - all data
test_that("toppFun save as xlsx works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = FALSE,
        format = "xlsx",
        verbose = FALSE
    )
    expect_true(file.exists(file.path(tmp_dir, "test_toppFun.xlsx")))
    df <- openxlsx::read.xlsx(file.path(tmp_dir, "test_toppFun.xlsx"), sheet = "toppData")
    expect_gt(nrow(df), 0)
})

## test save file as csv - split by cluster
test_that("toppFun save as split csv works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = TRUE,
        format = "csv",
        verbose = FALSE
    )
    expect_true(length(list.files(tmp_dir, pattern = "\\.csv$")) > 0)
    df <- read.csv(file.path(tmp_dir, "test_toppFun_CD8_T.csv"), header = TRUE)
    expect_gt(nrow(df), 0)
})

## test save file as csv - all data
test_that("toppFun save as csv works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = FALSE,
        format = "csv",
        verbose = FALSE
    )
    expect_true(file.exists(file.path(tmp_dir, "test_toppFun.csv")))
    df <- read.csv(file.path(tmp_dir, "test_toppFun.csv"), header = TRUE)
    expect_gt(nrow(df), 0)
})

## test save file as tsv - split by cluster
test_that("toppFun save as split tsv works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = TRUE,
        format = "tsv",
        verbose = FALSE
    )
    expect_true(length(list.files(tmp_dir, pattern = "\\.tsv$")) > 0)
    df <- read.table(file.path(tmp_dir, "test_toppFun_CD8_T.tsv"), header = TRUE, sep = "\t")
    expect_gt(nrow(df), 0)
})

## test save file as tsv - all data
test_that("toppFun save as tsv works", {
    data("ifnb.de")
    toppData <- toppFun(ifnb.de,
        type = "degs",
        topp_categories = "GeneOntologyMolecularFunction",
        cluster_col = "celltype",
        gene_col = "gene",
        logFC_col = "avg_log2FC",
        p_val_col = "p_val_adj",
        verbose = FALSE
    )

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppData,
        filename = "test_toppFun",
        save_dir = tmp_dir,
        split = FALSE,
        format = "tsv",
        verbose = FALSE
    )
    expect_true(file.exists(file.path(tmp_dir, "test_toppFun.tsv")))
    df <- read.table(file.path(tmp_dir, "test_toppFun.tsv"), header = TRUE, sep = "\t")
    expect_gt(nrow(df), 0)
})
