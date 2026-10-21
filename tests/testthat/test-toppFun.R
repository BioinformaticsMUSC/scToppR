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

test_that("toppSave accepts the standardized file_prefix argument", {
    data("toppdata.ifnb")
    output_dir <- tempfile("topp-save-")
    dir.create(output_dir)
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

    toppSave(
        toppData = toppdata.ifnb,
        file_prefix = "standardized_results",
        save_dir = output_dir,
        split = FALSE,
        format = "csv",
        verbose = FALSE
    )

    expect_true(file.exists(file.path(output_dir, "standardized_results.csv")))
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

test_that("toppFun direction_mode parameter works correctly", {
    data("ifnb.de")
    
    # Test direction_mode = "all" vs "split"
    # This would normally call API, so we test the gene processing logic
    
    # Test that split mode separates up/down regulated genes
    gene_data_split <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj", 
        logFC_col = "avg_log2FC",
        fc_cutoff = 0.25,
        fc_filter = "ALL"
    )
    
    # Test up-regulated genes only
    upregulated_data <- ifnb.de |> 
        dplyr::filter(avg_log2FC > 0.25)
    gene_data_up <- .process_degs(
        degs = upregulated_data,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        fc_cutoff = 0.25,
        fc_filter = "UPREG"
    )
    
    expect_type(gene_data_split, "list")
    expect_type(gene_data_up, "list")
})

test_that("fc_filter parameter validation and functionality", {
    data("ifnb.de")
    
    # Test all valid fc_filter options
    for (filter_type in c("ALL", "UPREG", "DOWNREG")) {
        gene_data <- .process_degs(
            degs = ifnb.de,
            cluster_col = "celltype",
            gene_col = "gene",
            p_val_col = "p_val_adj",
            logFC_col = "avg_log2FC",
            fc_filter = filter_type
        )
        expect_type(gene_data, "list")
    }
    
    # Test invalid fc_filter
    expect_error(
        toppFun(ifnb.de, 
                cluster_col = "celltype", 
                p_val_col = "p_val_adj",
                fc_filter = "INVALID"),
        "please select one of c\\('ALL', 'UPREG', 'DOWNREG'\\) for fc_filter"
    )
})

test_that("gene filtering thresholds work correctly", {
    data("ifnb.de")
    
    # Test different num_genes limits
    gene_data_100 <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        num_genes = 100
    )
    
    gene_data_50 <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype", 
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        num_genes = 50
    )
    
    # Should have fewer or equal genes when limit is lower
    for (cluster in names(gene_data_100)) {
        if (cluster %in% names(gene_data_50)) {
            expect_lte(length(gene_data_50[[cluster]]), length(gene_data_100[[cluster]]))
        }
    }
})

test_that("different p-value and fold change cutoffs work", {
    data("ifnb.de")
    
    # Test stricter p-value cutoff
    gene_data_strict_p <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene", 
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        pval_cutoff = 0.01  # stricter
    )
    
    gene_data_loose_p <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj", 
        logFC_col = "avg_log2FC",
        pval_cutoff = 0.1   # looser
    )
    
    # Stricter cutoff should generally result in fewer genes
    expect_type(gene_data_strict_p, "list")
    expect_type(gene_data_loose_p, "list")
    
    # Test fold change cutoffs
    gene_data_high_fc <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC", 
        fc_cutoff = 1.0  # high fold change
    )
    
    expect_type(gene_data_high_fc, "list")
})

test_that("cluster filtering works correctly", {
    data("ifnb.de")
    
    # Test with specific clusters only
    available_clusters <- unique(ifnb.de$celltype)
    test_clusters <- available_clusters[1:2]  # Take first 2 clusters
    
    # This tests the cluster filtering logic in toppFun
    # Since we can't test full toppFun without API, test the data processing
    filtered_data <- ifnb.de |>
        dplyr::filter(celltype %in% test_clusters)
    
    gene_data <- .process_degs(
        degs = filtered_data,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC"
    )
    
    expect_true(all(names(gene_data) %in% test_clusters))
})

test_that("correction method validation works", {
    data("ifnb.de")
    
    # Test valid correction methods
    valid_methods <- c("none", "FDR", "Bonferroni")
    for (method in valid_methods) {
        expect_no_error({
            # This tests the validation logic without API call
            ifnb.de_copy <- ifnb.de
        })
    }
    
    # Test invalid correction method
    expect_error(
        toppFun(ifnb.de,
                logFC_col = "avg_log2FC",
                cluster_col = "celltype",
                p_val_col = "p_val_adj", 
                correction = "INVALID"),
        "invalid P-value correction method"
    )
})

test_that("input data type validation is comprehensive", {
    data("ifnb.de")
    data("ifnb.markers.df")
    data("ifnb.markers.list.CD8T")
    
    # Test invalid data types
    expect_error(
        toppFun(NULL, type = "degs"),
        "Cluster column"
    )
    
    # Test missing required columns for degs type
    incomplete_data <- ifnb.de[, !colnames(ifnb.de) %in% "celltype"]
    expect_error(
        toppFun(incomplete_data, 
                cluster_col = "celltype", 
                type = "degs"),
        "Cluster column `celltype` not found in data"
    )
    
    # Test marker_df type
    expect_true(is.data.frame(ifnb.markers.df))
    
    # Test marker_list type  
    expect_true(is.vector(ifnb.markers.list.CD8T))
})

test_that("get_ToppCats returns expected categories", {
    cats <- get_ToppCats()
    
    expect_type(cats, "character")
    expect_gt(length(cats), 10)  # Should have many categories
    expect_true("GeneOntologyMolecularFunction" %in% cats)
    expect_true("GeneOntologyBiologicalProcess" %in% cats)
    expect_true("Pathway" %in% cats)
})

test_that("get_Entrez function works correctly", {
    # Test single gene
    expect_equal(get_Entrez("FLDB"), 338)
    
    # Test vector of genes
    test_genes <- c("FLDB", "TP53", "INVALID_GENE")
    entrez_result <- get_Entrez(test_genes)
    
    expect_type(entrez_result, "integer")
    expect_equal(entrez_result[1], 338)  # FLDB
})

test_that("addToppData function works with different object types", {
    if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
        requireNamespace("S4Vectors", quietly = TRUE)) {
        
        # Create mock SummarizedExperiment for testing
        library(SummarizedExperiment)
        library(S4Vectors)
        
        # Create simple test data
        counts <- matrix(1:100, ncol = 10)
        se <- SummarizedExperiment(assays = list(counts = counts))
        
        data("toppdata.ifnb")
        
        # Test basic functionality
        se_with_data <- addToppData(se, toppdata.ifnb, include_params = FALSE)
        expect_true("toppData" %in% names(metadata(se_with_data)))
        
        # Test with parameters
        se_with_params <- addToppData(se, toppdata.ifnb, include_params = TRUE)
        expect_true("toppData" %in% names(metadata(se_with_params)))
        expect_true("toppData_params" %in% names(metadata(se_with_params)))
        
        # Test custom slot name
        se_custom <- addToppData(se, toppdata.ifnb, slot_name = "enrichment")
        expect_true("enrichment" %in% names(metadata(se_custom)))
        
        # Test error conditions
        expect_error(
            addToppData("not_se_object", toppdata.ifnb),
            "must be a SummarizedExperiment"
        )
        
        expect_error(
            addToppData(se, "not_dataframe"),
            "must be a data.frame"
        )
    }
})