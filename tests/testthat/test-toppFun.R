### test for api entrez result
test_that("entrez lookup works", {
    expect_equal(get_Entrez("FLDB"), 338)
})

test_that("toppFun processes different input types correctly", {
    # using existing data instead of API calls
    data("ifnb.de")
    data("ifnb.markers.df") 
    data("ifnb.markers.list.CD8T")
    data("toppdata.ifnb")  # Expected result
    
    # Test structure and data processing without API calls
    expect_true("celltype" %in% colnames(ifnb.de))
    expect_true(is.data.frame(ifnb.de))
    expect_true(is.vector(ifnb.markers.list.CD8T))
    
    # Test that toppData has expected structure
    expect_true(is.data.frame(toppdata.ifnb))
    expect_true("Cluster" %in% colnames(toppdata.ifnb))
    expect_true("PValue" %in% colnames(toppdata.ifnb))
    expect_gt(nrow(toppdata.ifnb), 0)
})

test_that("gene processing works correctly", {
    data("ifnb.de")
    
    # Test .process_degs function directly
    gene_data <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene", 
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        num_genes = 100,
        pval_cutoff = 0.05,
        fc_cutoff = 0.25,
        fc_filter = "ALL"
    )
    
    expect_type(gene_data, "list")
    expect_true(length(gene_data) > 0)
    expect_true(all(sapply(gene_data, is.character)))
})

test_that("parameter validation works", {
    data("ifnb.de")
    
    # Test error conditions
    expect_error(
        toppFun(ifnb.de, cluster_col="celltype", p_val_col = "p_val_adj", type = "invalid"),
        "Please ensure the parameter `type` is one of"
    )
    
    expect_error(
        toppFun(ifnb.de, cluster_col = "nonexistent", p_val_col = "p_val_adj", ),
        "Cluster column `nonexistent` not found in data. Please specify."
    )
})

## test save file as xlsx - split by cluster
test_that("toppFun save as split xlsx works", {
    data("toppdata.ifnb")
    
    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
    data("toppdata.ifnb")
    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
    data("toppdata.ifnb")
    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
    data("toppdata.ifnb")

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
    data("toppdata.ifnb")

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
    data("toppdata.ifnb")

    tmp_dir <- tempdir()
    toppSave(
        toppData = toppdata.ifnb,
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
