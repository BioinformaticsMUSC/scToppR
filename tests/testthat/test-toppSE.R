skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")

make_degs <- function() {
    data("ifnb.de", envir = environment())
    ifnb.de[ifnb.de$celltype %in% c("CD14 Mono", "CD8 T", "NK"), ]
}

make_degs_se <- function(degs = make_degs()) {
    toppSE(degs,
        cluster_col = "celltype", gene_col = "gene",
        p_val_col = "p_val_adj", logFC_col = "avg_log2FC"
    )
}

# stand-in for the ToppGene query: one fake hit per cluster (and direction)
mock_results <- function(directions = NULL) {
    data("toppdata.ifnb", envir = environment())
    res <- toppdata.ifnb[toppdata.ifnb$Cluster %in% c("CD14 Mono", "CD8 T"), ]
    if (!is.null(directions)) {
        res$Direction <- ifelse(seq_len(nrow(res)) %% 2 == 0, "up", "down")
    }
    res
}

test_that("toppSE builds a container from DEGs with standardized columns", {
    degs <- make_degs()
    se <- make_degs_se(degs)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(nrow(se), 0)
    expect_setequal(colnames(se), c("CD14 Mono", "CD8 T", "NK"))

    stored <- S4Vectors::metadata(se)$DEGs
    expect_named(stored, c("cluster", "gene", "logFC", "p_val"))
    expect_equal(nrow(stored), nrow(degs))

    col_data <- SummarizedExperiment::colData(se)
    expect_equal(as.integer(col_data["NK", "n_input_genes"]), sum(degs$celltype == "NK"))
    expect_true(all(col_data$Direction == "all"))
})

test_that("toppSE validates its input", {
    degs <- make_degs()
    expect_error(toppSE(degs), "not found in the DEG data")
    expect_error(
        toppSE(degs, cluster_col = "celltype", type = "bogus"),
        "`type` is one of"
    )
    expect_error(
        toppSE(degs[0, ], cluster_col = "celltype"),
        "No genes found"
    )
})

test_that("toppSE accepts marker lists and marker data.frames", {
    markers <- list(A = c("IFNG", "FOXP3"), B = c("CD8A", NA, ""))
    se <- toppSE(markers)
    stored <- S4Vectors::metadata(se)$DEGs
    expect_equal(colnames(se), c("A", "B"))
    expect_equal(stored$gene, c("IFNG", "FOXP3", "CD8A"))
    expect_true(all(is.na(stored$logFC)))
    expect_equal(S4Vectors::metadata(se)$toppInput$type, "marker_df")

    se_vec <- toppSE(c("IFNG", "FOXP3"))
    expect_equal(S4Vectors::metadata(se_vec)$toppInput$type, "marker_list")
    expect_equal(colnames(se_vec), "genes")
})

test_that("toppSE converts toppFun results and toppData() round-trips them", {
    data("toppdata.ifnb")
    se <- toppSE(toppdata.ifnb)

    expect_equal(S4Vectors::metadata(se)$toppInput$type, "results")
    expect_setequal(colnames(se), unique(toppdata.ifnb$Cluster))
    expect_setequal(
        SummarizedExperiment::assayNames(se),
        c("PValue", "QValueFDRBH", "QValueFDRBY", "QValueBonferroni", "GenesInQuery", "GenesInTermInQuery", "Genes")
    )

    back <- toppData(se)
    sort_df <- function(x) {
        x <- x[order(x$Cluster, x$Category, x$ID, x$Name), colnames(toppdata.ifnb)]
        rownames(x) <- NULL
        x
    }
    expect_equal(sort_df(back), sort_df(toppdata.ifnb))
})

test_that("toppSE results input keeps a Direction column for split results", {
    res <- mock_results(directions = TRUE)
    se <- toppSE(res)
    expect_true(all(c("CD14 Mono_up", "CD14 Mono_down") %in% colnames(se)))
    expect_true("Direction" %in% colnames(toppData(se)))
    expect_equal(nrow(toppData(se)), nrow(res))
})

test_that("toppData reads legacy metadata slots and validates input", {
    data("toppdata.ifnb")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, 2, 2))
    )
    se <- addToppData(se, toppdata.ifnb, include_params = FALSE)
    expect_identical(toppData(se), toppdata.ifnb)
    expect_error(toppData(se, slot_name = "nope"), "not found in metadata slot 'nope'")
    expect_error(toppData("x"), "must be a SummarizedExperiment")
})

test_that("addToppData records the clusters of the results", {
    data("toppdata.ifnb")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, 2, 2))
    )
    se <- addToppData(se, toppdata.ifnb)
    params <- S4Vectors::metadata(se)$toppData_params
    expect_setequal(params$clusters, unique(toppdata.ifnb$Cluster))
})

test_that("toppFun on a toppSE uses the stored columns and fills the assays", {
    se <- make_degs_se()
    captured <- NULL
    testthat::local_mocked_bindings(
        toppFun.default = function(input_data, ...) {
            captured <<- list(input_data = input_data, args = list(...))
            mock_results()
        },
        .package = "scToppR"
    )

    out <- toppFun(se, num_genes = 50, fc_cutoff = 0.5)

    expect_equal(captured$args$cluster_col, "cluster")
    expect_equal(captured$args$gene_col, "gene")
    expect_equal(captured$args$logFC_col, "logFC")
    expect_equal(captured$args$p_val_col, "p_val")
    expect_equal(captured$args$num_genes, 50)
    expect_equal(captured$args$key_type, "SYMBOL")
    expect_named(captured$input_data, c("cluster", "gene", "logFC", "p_val"))

    expect_gt(nrow(out), 0)
    expect_setequal(colnames(out), c("CD14 Mono", "CD8 T", "NK"))
    # NK had no results but keeps its column
    expect_true(all(is.na(SummarizedExperiment::assay(out, "PValue")[, "NK"])))
    params <- S4Vectors::metadata(out)$toppParams
    expect_equal(params$num_genes, 50)
    expect_equal(params$fc_cutoff, 0.5)
    expect_equal(nrow(toppData(out)), nrow(mock_results()))
    # the DEGs stay in the object
    expect_equal(nrow(S4Vectors::metadata(out)$DEGs), nrow(make_degs()))
})

test_that("toppFun on a toppSE supports clusters, split mode and re-running", {
    se <- make_degs_se()
    testthat::local_mocked_bindings(
        toppFun.default = function(input_data, direction_mode, ...) {
            if (direction_mode == "split") mock_results(directions = TRUE) else mock_results()
        },
        .package = "scToppR"
    )

    split_se <- toppFun(se, direction_mode = "split", clusters = c("CD14 Mono", "CD8 T"))
    expect_setequal(
        colnames(split_se),
        c("CD14 Mono_up", "CD14 Mono_down", "CD8 T_up", "CD8 T_down")
    )
    col_data <- SummarizedExperiment::colData(split_se)
    expect_setequal(col_data$Direction, c("up", "down"))

    expect_error(toppFun(se, clusters = "nope"), "not found in the stored DEGs")

    first <- toppFun(se)
    expect_warning(toppFun(first), "Replacing the enrichment results")
})

test_that("toppFun on a toppSE validates the object and options", {
    se <- make_degs_se()
    expect_error(toppFun(se, direction_mode = "bogus"), "direction_mode")

    marker_se <- toppSE(list(A = c("IFNG", "FOXP3")))
    expect_error(toppFun(marker_se, direction_mode = "split"), "requires DEGs")

    plain <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1, 2, 2)))
    expect_error(toppFun(plain), "created with toppSE")

    results_se <- toppSE(mock_results())
    expect_error(toppFun(results_se), "Add them with addDEGs")
})

test_that("toppFun on marker input queries per-cluster gene lists", {
    se <- toppSE(list(A = c("IFNG", "FOXP3"), B = c("CD8A", "CD4")))
    captured <- NULL
    testthat::local_mocked_bindings(
        toppFun.default = function(input_data, type, ...) {
            captured <<- list(input_data = input_data, type = type)
            res <- mock_results()
            res$Cluster <- ifelse(res$Cluster == "CD14 Mono", "A", "B")
            res
        },
        .package = "scToppR"
    )
    toppFun(se)
    expect_equal(captured$type, "marker_df")
    expect_equal(captured$input_data, list(A = c("IFNG", "FOXP3"), B = c("CD8A", "CD4")))
})

test_that("toppFun dispatches for SingleCellExperiment-based containers", {
    skip_if_not_installed("SingleCellExperiment")
    se <- make_degs_se()
    sce <- methods::as(se, "SingleCellExperiment")
    expect_true(.is_toppSE(sce))
    testthat::local_mocked_bindings(
        toppFun.default = function(input_data, ...) mock_results(),
        .package = "scToppR"
    )
    out <- toppFun(sce)
    expect_gt(nrow(out), 0)
})

test_that("toppFun still dispatches to the default method for data.frames", {
    expect_true(is.function(getS3method("toppFun", "default")))
    expect_error(
        toppFun(make_degs(), type = "bogus"),
        "Please ensure the parameter `type`"
    )
})

test_that("addDEGs stores standardized DEGs and checks the clusters", {
    data("toppdata.ifnb")
    se <- toppSE(toppdata.ifnb)
    degs <- make_degs()

    expect_warning(
        out <- addDEGs(se, degs,
            cluster_col = "celltype", gene_col = "gene",
            p_val_col = "p_val_adj", logFC_col = "avg_log2FC"
        ),
        "No DEGs found for cluster"
    )
    stored <- S4Vectors::metadata(out)$DEGs
    expect_named(stored, c("cluster", "gene", "logFC", "p_val"))
    col_data <- SummarizedExperiment::colData(out)
    expect_equal(as.integer(col_data["NK", "n_input_genes"]), sum(degs$celltype == "NK"))

    renamed <- degs
    renamed$celltype <- paste0("x", renamed$celltype)
    expect_error(
        addDEGs(se, renamed, cluster_col = "celltype"),
        "None of the clusters"
    )
    expect_error(addDEGs(se, degs), "not found in the DEG data")
    expect_error(addDEGs("x", degs), "created with toppSE")
})

test_that("termGenes joins hit genes with DEG statistics", {
    data("toppdata.ifnb")
    data("ifnb.de")
    se <- toppSE(toppdata.ifnb)
    se <- suppressWarnings(addDEGs(se, ifnb.de,
        cluster_col = "celltype", gene_col = "gene",
        p_val_col = "p_val_adj", logFC_col = "avg_log2FC"
    ))

    out <- termGenes(se, "cytokine receptor binding", cluster = "CD14 Mono")
    expect_named(out, c("Cluster", "Category", "ID", "Name", "gene", "logFC", "p_val"))
    expected_genes <- strsplit(
        toppdata.ifnb$Genes[toppdata.ifnb$Name == "cytokine receptor binding" &
            toppdata.ifnb$Cluster == "CD14 Mono"][1], ",\\s*"
    )[[1]]
    expect_equal(out$gene, expected_genes)
    expect_true(any(!is.na(out$logFC)))

    # same statistics as the DEG table
    check <- ifnb.de[ifnb.de$celltype == "CD14 Mono" & ifnb.de$gene == out$gene[1], ]
    expect_equal(out$logFC[1], check$avg_log2FC[1])

    # lookup by ID, and NA statistics without DEGs
    no_degs <- termGenes(toppSE(toppdata.ifnb), "GO:0005126", cluster = "CD14 Mono")
    expect_true(all(is.na(no_degs$logFC)))

    expect_error(termGenes(se, "not a term"), "No enrichment results")
    expect_error(termGenes("x", "GO:0005126"), "created with toppSE")
})

test_that("plotting and saving work directly on a toppSE", {
    data("toppdata.pbmc")
    se <- toppSE(toppdata.pbmc)
    cluster <- unique(toppdata.pbmc$Cluster)[1]

    plot <- suppressWarnings(toppPlot(
        se,
        category = "GeneOntologyMolecularFunction", clusters = cluster, save = FALSE
    ))
    expect_s3_class(plot, "ggplot")

    balloon <- toppBalloon(se, categories = "GeneOntologyMolecularFunction", save = FALSE)
    expect_s3_class(balloon, "ggplot")

    output_dir <- tempfile("topp-se-")
    dir.create(output_dir)
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
    toppSave(se, file_prefix = "se_results", save_dir = output_dir, split = FALSE, format = "csv", verbose = FALSE)
    saved <- read.csv(file.path(output_dir, "se_results.csv"))
    expect_equal(nrow(saved), nrow(toppdata.pbmc))
})

test_that("toppVolc highlights term genes from a toppSE and from a data.frame", {
    skip_if_not_installed("ggrepel")
    data("toppdata.ifnb")
    data("ifnb.de")
    cat_name <- "GeneOntologyMolecularFunction"

    se <- toppSE(toppdata.ifnb)
    se <- suppressWarnings(addDEGs(se, ifnb.de,
        cluster_col = "celltype", gene_col = "gene",
        p_val_col = "p_val_adj", logFC_col = "avg_log2FC"
    ))
    plots <- toppVolc(se, category = cat_name, cluster = "NK", num_terms = 2)
    expect_length(plots, 2)
    expect_s3_class(plots[[1]], "ggplot")
    plot_data <- plots[[1]]$data
    expect_equal(nrow(plot_data), sum(ifnb.de$celltype == "NK" & !is.na(ifnb.de$p_val_adj)))
    expect_true(any(plot_data$Highlight == "In gene set"))

    one <- toppVolc(toppdata.ifnb,
        degs = ifnb.de, category = cat_name, cluster = "NK", num_terms = 1,
        deg_cluster_col = "celltype", p_val_col = "p_val_adj"
    )
    expect_s3_class(one, "ggplot")
    expect_s3_class(
        toppVolc(se, category = cat_name, cluster = "NK", num_terms = 2, combine = TRUE),
        "patchwork"
    )

    out_dir <- tempfile("volc-")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    toppVolc(se, category = cat_name, cluster = "NK", num_terms = 1, save = TRUE, save_dir = out_dir)
    expect_length(list.files(out_dir, pattern = "^toppVolc_NK_.*\\.pdf$"), 1)

    expect_error(toppVolc(se, cluster = "NK"), "`category`")
    expect_error(toppVolc(se, category = cat_name), "exactly one `cluster`")
    expect_error(toppVolc(se, category = cat_name, cluster = "nope"), "not found in the results")
    expect_error(toppVolc(toppdata.ifnb, category = cat_name, cluster = "NK"), "`degs` is required")
})
