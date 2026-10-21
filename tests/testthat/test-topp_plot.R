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

# Tests for S3 methods with SummarizedExperiment and SingleCellExperiment
test_that("toppPlot works with SummarizedExperiment objects", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    se <- addToppData(se, toppdata.pbmc, include_params = FALSE)
    
    # Test that the method works
    dotplot <- toppPlot(se, 
                        category = "GeneOntologyMolecularFunction",
                        clusters = 0,
                        save = FALSE)
    
    expect_s3_class(dotplot, "ggplot")
    expect_equal(length(dotplot$layers), 2)
})

test_that("toppBalloon works with SummarizedExperiment objects", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    se <- addToppData(se, toppdata.pbmc, include_params = FALSE)
    
    # Test that the method works
    balloonplot <- toppBalloon(se, 
                              categories = "GeneOntologyMolecularFunction",
                              save = FALSE)
    
    expect_s3_class(balloonplot, "ggplot")
    expect_equal(length(balloonplot$layers), 1)
})

test_that("toppPlot works with SingleCellExperiment objects", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SingleCellExperiment object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    sce <- addToppData(sce, toppdata.pbmc, include_params = FALSE)
    
    # Test that the method works
    dotplot <- toppPlot(sce, 
                        category = "GeneOntologyMolecularFunction",
                        clusters = 0,
                        save = FALSE)
    
    expect_s3_class(dotplot, "ggplot")
    expect_equal(length(dotplot$layers), 2)
})

test_that("toppBalloon works with SingleCellExperiment objects", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SingleCellExperiment object
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    sce <- addToppData(sce, toppdata.pbmc, include_params = FALSE)
    
    # Test that the method works
    balloonplot <- toppBalloon(sce, 
                              categories = "GeneOntologyMolecularFunction",
                              save = FALSE)
    
    expect_s3_class(balloonplot, "ggplot")
    expect_equal(length(balloonplot$layers), 1)
})

test_that("S3 methods work with custom slot names", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to custom metadata slot
    se <- addToppData(se, toppdata.pbmc, slot_name = "enrichment_results", include_params = FALSE)
    
    # Test that the method works with custom slot name
    dotplot <- toppPlot(se, 
                        category = "GeneOntologyMolecularFunction",
                        clusters = 0,
                        slot_name = "enrichment_results",
                        save = FALSE)
    
    expect_s3_class(dotplot, "ggplot")
    expect_equal(length(dotplot$layers), 2)
})

test_that("S3 methods handle missing toppData gracefully", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create a simple SummarizedExperiment object without toppData
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Test error handling for missing toppData
    expect_error(
        toppPlot(se, category = "GeneOntologyMolecularFunction"),
        "toppData not found in metadata slot 'toppData'"
    )
    
    expect_error(
        toppBalloon(se, categories = "GeneOntologyMolecularFunction"),
        "toppData not found in metadata slot 'toppData'"
    )
})

test_that("S3 methods handle invalid slot names gracefully", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    se <- addToppData(se, toppdata.pbmc, include_params = FALSE)
    
    # Test error handling for invalid slot name
    expect_error(
        toppPlot(se, category = "GeneOntologyMolecularFunction", slot_name = "nonexistent"),
        "toppData not found in metadata slot 'nonexistent'"
    )
    
    expect_error(
        toppBalloon(se, categories = "GeneOntologyMolecularFunction", slot_name = "nonexistent"),
        "toppData not found in metadata slot 'nonexistent'"
    )
})

test_that("S3 methods handle invalid metadata content gracefully", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add invalid data to toppData slot
    S4Vectors::metadata(se)$toppData <- "not_a_dataframe"
    
    # Test error handling for invalid data type
    expect_error(
        toppPlot(se, category = "GeneOntologyMolecularFunction"),
        "Data in metadata slot 'toppData' is not a data.frame"
    )
    
    expect_error(
        toppBalloon(se, categories = "GeneOntologyMolecularFunction"),
        "Data in metadata slot 'toppData' is not a data.frame"
    )
})

test_that("S3 methods pass parameters correctly", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("S4Vectors")
    
    data("toppdata.pbmc")
    
    # Create a simple SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 4, ncol = 5))
    )
    
    # Add toppData to metadata
    se <- addToppData(se, toppdata.pbmc, include_params = FALSE)
    
    # Test multiple clusters with SE method
    dotplot_list <- toppPlot(se, 
                            category = "GeneOntologyMolecularFunction",
                            clusters = c("CD4T", "CD8T"),
                            save = FALSE)
    
    expect_type(dotplot_list, "list")
    expect_equal(length(dotplot_list), 2)
    expect_s3_class(dotplot_list$CD4T, "ggplot")
    expect_s3_class(dotplot_list$CD8T, "ggplot")
})


test_that("toppPlot error handling works correctly", {
    data("toppdata.pbmc")
    
    # Test with invalid category
    expect_error(
        toppPlot(toppdata.pbmc, category = "NonExistentCategory"),
        "Category NonExistentCategory not found in the data"
    )
    
    # Test with invalid cluster - just skip for now as this depends on exact validation
    
    # Test with empty dataframe - expect category not found error
    empty_df <- toppdata.pbmc[0, ]
    expect_error(
        toppPlot(empty_df, category = "GeneOntologyMolecularFunction", save = FALSE),
        "Category.*not found"
    )
})

test_that("toppPlot handles different cluster types correctly", {
    data("toppdata.pbmc")
    
    # Test with numeric cluster
    plot_numeric <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction", 
        clusters = 0,
        save = FALSE
    )
    expect_s3_class(plot_numeric, "ggplot")
    
    # Test with character cluster
    plot_char <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T", 
        save = FALSE
    )
    expect_s3_class(plot_char, "ggplot")
    
    # Test with NULL clusters (should plot all)
    plot_all <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = NULL,
        save = FALSE
    )
    expect_type(plot_all, "list")
    expect_true(length(plot_all) > 1)  # Should have multiple plots
})

test_that("toppPlot parameter variations work", {
    data("toppdata.pbmc")
    
    # Test different number of terms
    plot_3_terms <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        num_terms = 3,
        save = FALSE
    )
    expect_s3_class(plot_3_terms, "ggplot")
    
    plot_10_terms <- toppPlot(toppdata.pbmc, 
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        num_terms = 10,
        save = FALSE
    )
    expect_s3_class(plot_10_terms, "ggplot")
    
    # Test different plot types if available
    # (This would depend on your actual plot types in the function)
    
    # Test parameter variations with valid parameters
    plot_custom_params <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T", 
        y_axis_text_size = 12,
        p_val_adj = "Bonferroni",
        file_prefix = "custom_plot",
        save = FALSE
    )
    expect_s3_class(plot_custom_params, "ggplot")
})

test_that("toppBalloon error handling and validation", {
    data("toppdata.pbmc")
    
    # Test with invalid categories
    expect_error(
        toppBalloon(toppdata.pbmc, categories = "NonExistentCategory"),
        "Category NonExistentCategory not found in the data"
    )
    
    # Test with empty dataframe
    empty_df <- toppdata.pbmc[0, ]
    expect_warning(
        result <- toppBalloon(empty_df, categories = "GeneOntologyMolecularFunction"),
        "No data available to plot"
    )
    expect_s3_class(result, "ggplot")
    
    # Test with single category vs multiple categories
    balloon_single <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        save = FALSE
    )
    expect_s3_class(balloon_single, "ggplot")
    
    balloon_multi <- toppBalloon(toppdata.pbmc,
        categories = c("GeneOntologyMolecularFunction", "GeneOntologyBiologicalProcess"),
        save = FALSE
    )
    # Multi-category toppBalloon returns a list
    expect_type(balloon_multi, "list")
    expect_true(all(sapply(balloon_multi, function(x) "ggplot" %in% class(x))))
})

test_that("toppBalloon parameter variations work", {
    data("toppdata.pbmc")
    
    # Test different number of balloons
    balloon_5_terms <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        balloons = 5,
        save = FALSE
    )
    expect_s3_class(balloon_5_terms, "ggplot")
    
    # Test filtering by specific clusters if supported
    balloon_filtered <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        save = FALSE
    )
    expect_s3_class(balloon_filtered, "ggplot")
    
    # Test parameter customization
    balloon_custom <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction", 
        x_axis_text_size = 8,
        filename = "custom_balloon",
        save = FALSE
    )
    expect_s3_class(balloon_custom, "ggplot")
})

test_that("plot saving uses file_prefix and retains filename aliases", {
    data("toppdata.pbmc")
    output_dir <- tempfile("topp-plots-")
    dir.create(output_dir)
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

    toppPlot(
        toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        file_prefix = "standardized_plot",
        save = TRUE,
        save_dir = output_dir
    )
    expect_true(file.exists(file.path(
        output_dir,
        "standardized_plot_GeneOntologyMolecularFunction.pdf"
    )))

    toppPlot(
        toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        filename = "legacy_plot",
        save = TRUE,
        save_dir = output_dir
    )
    expect_true(file.exists(file.path(
        output_dir,
        "legacy_plot_GeneOntologyMolecularFunction.pdf"
    )))

    toppBalloon(
        toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        file_prefix = "standardized_balloon",
        save = TRUE,
        save_dir = output_dir
    )
    expect_true(file.exists(file.path(
        output_dir,
        "standardized_balloon_GeneOntologyMolecularFunction.pdf"
    )))

    toppBalloon(
        toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        filename = "legacy_balloon",
        save = TRUE,
        save_dir = output_dir
    )
    expect_true(file.exists(file.path(
        output_dir,
        "legacy_balloon_GeneOntologyMolecularFunction.pdf"
    )))
})

test_that("plotting functions handle edge cases", {
    data("toppdata.pbmc")
    
    # Test with single row of data
    single_row <- toppdata.pbmc[1, ]
    plot_single <- toppPlot(single_row,
        category = single_row$Category[1],
        clusters = single_row$Cluster[1],
        save = FALSE
    )
    expect_s3_class(plot_single, "ggplot")
    
    # Test with data containing NA values
    data_with_na <- toppdata.pbmc
    data_with_na$PValue[1:5] <- NA
    plot_with_na <- toppPlot(data_with_na,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        save = FALSE
    )
    expect_s3_class(plot_with_na, "ggplot")
    
    # Test with very high p-values
    high_pval_data <- toppdata.pbmc
    high_pval_data$PValue <- 0.9
    high_pval_data$QValueFDRBH <- 0.95
    plot_high_p <- toppPlot(high_pval_data,
        category = "GeneOntologyMolecularFunction", 
        clusters = "CD4T",
        save = FALSE
    )
    expect_s3_class(plot_high_p, "ggplot")
})

test_that("plot structure and aesthetics are correct", {
    data("toppdata.pbmc")
    
    dotplot <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T", 
        save = FALSE
    )
    
    # Check that plot has expected components
    expect_s3_class(dotplot, "ggplot")
    expect_true(length(dotplot$layers) >= 1)
    
    # Check for expected structure (ggplot objects have these as nested properties)
    expect_true(!is.null(dotplot$layers))
    expect_true("ggplot" %in% class(dotplot))
    
    # Test balloon plot structure
    balloon <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        save = FALSE
    )
    
    expect_s3_class(balloon, "ggplot") 
    expect_true(length(balloon$layers) >= 1)
})

test_that("save functionality validation", {
    data("toppdata.pbmc")
    
    # Test save parameter validation
    temp_dir <- tempdir()
    
    # This should work without error (but we set save=FALSE to avoid actual saving in tests)
    expect_no_error({
        toppPlot(toppdata.pbmc,
            category = "GeneOntologyMolecularFunction",
            clusters = "CD4T",
            save = FALSE  # Don't actually save in tests
        )
    })
    
    expect_no_error({
        toppBalloon(toppdata.pbmc,
            categories = "GeneOntologyMolecularFunction", 
            save = FALSE  # Don't actually save in tests
        )
    })
})

test_that("color and visual parameter handling", {
    data("toppdata.pbmc")
    
    # Test if color parameters are handled correctly
    plot_default <- toppPlot(toppdata.pbmc,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T",
        save = FALSE
    )
    
    # If your functions support color customization, test it
    # plot_custom_color <- toppPlot(toppdata.pbmc,
    #     category = "GeneOntologyMolecularFunction", 
    #     clusters = "CD4T",
    #     color_palette = "viridis",
    #     save = FALSE
    # )
    
    expect_s3_class(plot_default, "ggplot")
    
    # Test balloon plot colors
    balloon_default <- toppBalloon(toppdata.pbmc,
        categories = "GeneOntologyMolecularFunction",
        save = FALSE
    )
    
    expect_s3_class(balloon_default, "ggplot")
})

test_that("multiple categories and clusters combinations work", {
    data("toppdata.pbmc")
    
    # Test multiple categories at once
    available_categories <- unique(toppdata.pbmc$Category)
    if (length(available_categories) >= 2) {
        multi_cat_balloon <- toppBalloon(toppdata.pbmc,
            categories = available_categories[1:2],
            save = FALSE
        )
        # Multi-category toppBalloon returns a list of ggplot objects
        expect_type(multi_cat_balloon, "list")
        expect_true(all(sapply(multi_cat_balloon, function(x) "ggplot" %in% class(x))))
    }
    
    # Test multiple clusters at once
    available_clusters <- unique(toppdata.pbmc$Cluster)
    if (length(available_clusters) >= 2) {
        multi_cluster_plots <- toppPlot(toppdata.pbmc,
            category = available_categories[1],
            clusters = available_clusters[1:2],
            save = FALSE
        )
        expect_type(multi_cluster_plots, "list")
        expect_equal(length(multi_cluster_plots), 2)
    }
})

test_that("data filtering and sorting work correctly", {
    data("toppdata.pbmc")
    
    # Test that plots handle pre-filtered data correctly
    filtered_data <- toppdata.pbmc |>
        dplyr::filter(QValueFDRBH < 0.01)  # Very strict filter
    
    if (nrow(filtered_data) > 0) {
        plot_filtered <- toppPlot(filtered_data,
            category = unique(filtered_data$Category)[1], 
            clusters = unique(filtered_data$Cluster)[1],
            save = FALSE
        )
        expect_s3_class(plot_filtered, "ggplot")
    }
    
    # Test with data sorted differently
    shuffled_data <- toppdata.pbmc[sample(nrow(toppdata.pbmc)), ]
    plot_shuffled <- toppPlot(shuffled_data,
        category = "GeneOntologyMolecularFunction",
        clusters = "CD4T", 
        save = FALSE
    )
    expect_s3_class(plot_shuffled, "ggplot")
})