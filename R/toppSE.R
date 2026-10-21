## Standard column sets used when converting between the long toppFun()
## data.frame and the term x cluster SummarizedExperiment representation
.topp_result_cols <- c(
    "Category", "ID", "Name", "PValue", "QValueFDRBH", "QValueFDRBY",
    "QValueBonferroni", "TotalGenes", "GenesInTerm", "GenesInQuery",
    "GenesInTermInQuery", "Source", "URL", "Genes"
)
.topp_assay_fields <- c(
    "PValue", "QValueFDRBH", "QValueFDRBY", "QValueBonferroni",
    "GenesInQuery", "GenesInTermInQuery", "Genes"
)
.topp_term_fields <- c("Category", "ID", "Name", "TotalGenes", "GenesInTerm", "Source", "URL")
.topp_required_cols <- c("Category", "ID", "Name", "PValue", "GenesInTermInQuery")

.check_se_pkgs <- function() {
    for (pkg in c("SummarizedExperiment", "S4Vectors")) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(pkg, " package is required but not available. ",
                "Please install it with: BiocManager::install('", pkg, "')",
                call. = FALSE
            )
        }
    }
}

.is_toppSE <- function(x) {
    inherits(x, "SummarizedExperiment") && isTRUE(S4Vectors::metadata(x)$toppSE)
}

.looks_like_results <- function(x) {
    is.data.frame(x) && all(.topp_required_cols %in% colnames(x))
}

.topp_col_names <- function(cluster, direction) {
    as.character(ifelse(direction == "all", cluster, paste(cluster, direction, sep = "_")))
}

.empty_results <- function() {
    chr <- c("Category", "ID", "Name", "Source", "URL", "Genes")
    out <- lapply(.topp_result_cols, function(col) if (col %in% chr) character() else numeric())
    names(out) <- .topp_result_cols
    as.data.frame(out, stringsAsFactors = FALSE)
}

## Standardize DEG results to the fixed columns cluster, gene, logFC, p_val
.standardize_degs <- function(degs, cluster_col, gene_col, logFC_col, p_val_col) {
    degs <- as.data.frame(degs)
    needed <- c(cluster_col, gene_col, logFC_col, p_val_col)
    missing_cols <- setdiff(needed, colnames(degs))
    if (length(missing_cols) > 0) {
        stop("Column(s) not found in the DEG data: ",
            paste0("`", missing_cols, "`", collapse = ", "),
            ". Please specify them with the `cluster_col`, `gene_col`, `logFC_col` and `p_val_col` arguments.",
            call. = FALSE
        )
    }
    data.frame(
        cluster = as.character(degs[[cluster_col]]),
        gene = as.character(degs[[gene_col]]),
        logFC = as.numeric(degs[[logFC_col]]),
        p_val = as.numeric(degs[[p_val_col]]),
        stringsAsFactors = FALSE
    )
}

## Marker lists carry no statistics, so logFC and p_val are NA
.standardize_markers <- function(markers, type) {
    if (type == "marker_list") {
        markers <- list(genes = markers)
    } else {
        markers <- as.list(markers)
        if (is.null(names(markers)) || any(names(markers) == "")) {
            stop("Marker input must have named columns/elements giving the cluster names.", call. = FALSE)
        }
    }
    long <- lapply(names(markers), function(cl) {
        genes <- as.character(markers[[cl]])
        genes <- genes[!is.na(genes) & nzchar(genes)]
        data.frame(
            cluster = rep(cl, length(genes)), gene = genes,
            logFC = rep(NA_real_, length(genes)), p_val = rep(NA_real_, length(genes)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, long)
}

.count_input_genes <- function(degs, clusters) {
    as.integer(table(factor(degs$cluster, levels = unique(clusters)))[clusters])
}

## Build the term x cluster SummarizedExperiment.
## `results` is the long toppFun() data.frame, `col_df` defines the columns
## (Cluster, Direction, n_input_genes) so clusters without hits are retained.
.build_topp_se <- function(results, col_df, md) {
    .check_se_pkgs()
    col_names <- .topp_col_names(col_df$Cluster, col_df$Direction)
    if (anyDuplicated(col_names)) {
        stop("Cluster/direction combinations must be unique.", call. = FALSE)
    }

    results <- if (is.null(results) || nrow(results) == 0) {
        .empty_results()
    } else {
        as.data.frame(results, stringsAsFactors = FALSE)
    }

    direction <- if ("Direction" %in% colnames(results)) as.character(results$Direction) else rep("all", nrow(results))
    col_key <- .topp_col_names(as.character(results$Cluster), direction)
    term_key <- paste(results$Category, results$ID, results$Name, sep = "|")
    if (!all(col_key %in% col_names)) {
        stop("Results contain clusters that are not defined in the object.", call. = FALSE)
    }

    dup <- duplicated(cbind(term_key, col_key))
    if (any(dup)) {
        warning(sum(dup), " duplicated term/cluster result(s) were dropped.", call. = FALSE)
        results <- results[!dup, , drop = FALSE]
        term_key <- term_key[!dup]
        col_key <- col_key[!dup]
    }

    keys <- unique(term_key)
    idx <- cbind(match(term_key, keys), match(col_key, col_names))

    assay_list <- list()
    for (field in intersect(.topp_assay_fields, colnames(results))) {
        blank <- if (field == "Genes") NA_character_ else NA_real_
        mat <- matrix(blank, nrow = length(keys), ncol = length(col_names), dimnames = list(keys, col_names))
        if (nrow(idx) > 0) {
            mat[idx] <- results[[field]]
        }
        assay_list[[field]] <- mat
    }

    row_df <- results[match(keys, term_key), intersect(.topp_term_fields, colnames(results)), drop = FALSE]
    rownames(row_df) <- keys

    col_df$n_query_genes <- NA_integer_
    if ("GenesInQuery" %in% colnames(results) && nrow(results) > 0) {
        q <- tapply(results$GenesInQuery, factor(col_key, levels = col_names), max)
        col_df$n_query_genes <- as.integer(q)
    }
    rownames(col_df) <- col_names

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = assay_list,
        rowData = S4Vectors::DataFrame(row_df, check.names = FALSE),
        colData = S4Vectors::DataFrame(col_df, check.names = FALSE)
    )
    md$toppSE <- TRUE
    S4Vectors::metadata(se) <- md
    se
}

## Convert a toppSE back to the long toppFun() data.frame
.se_to_results <- function(se) {
    p_val <- SummarizedExperiment::assay(se, "PValue")
    hit <- which(!is.na(p_val), arr.ind = TRUE)
    row_i <- hit[, 1]
    col_i <- hit[, 2]

    row_df <- as.data.frame(SummarizedExperiment::rowData(se))
    col_df <- as.data.frame(SummarizedExperiment::colData(se))
    assay_names <- SummarizedExperiment::assayNames(se)

    out <- row_df[row_i, intersect(.topp_term_fields, colnames(row_df)), drop = FALSE]
    rownames(out) <- NULL
    for (field in intersect(.topp_assay_fields, assay_names)) {
        out[[field]] <- SummarizedExperiment::assay(se, field)[hit]
    }
    out <- out[, intersect(.topp_result_cols, colnames(out)), drop = FALSE]
    out$Cluster <- col_df$Cluster[col_i]
    if (any(col_df$Direction != "all")) {
        out$Direction <- col_df$Direction[col_i]
    }

    category_order <- match(out$Category, unique(row_df$Category))
    out <- out[order(col_i, category_order, out$PValue), , drop = FALSE]
    rownames(out) <- NULL
    out
}

#' Create a term-by-cluster SummarizedExperiment for ToppGene analysis
#'
#' @description
#' `toppSE()` builds a SummarizedExperiment that holds everything for a scToppR
#' analysis in one place. Differential expression results (or marker genes) are
#' stored with standardized column names, so [toppFun()] can be run on the object
#' without specifying any column names. After [toppFun()] is run, the enrichment
#' results are stored with terms as rows and clusters as columns.
#'
#' @details
#' The object has the following structure:
#' \itemize{
#'   \item **Rows**: enriched terms, named `Category|ID|Name`.
#'   \item **Columns**: clusters (or `cluster_up`/`cluster_down` when
#'     `direction_mode = "split"`). `colData` has the columns `Cluster`,
#'     `Direction`, `n_input_genes` and `n_query_genes`.
#'   \item **rowData**: term-level information (`Category`, `ID`, `Name`,
#'     `TotalGenes`, `GenesInTerm`, `Source`, `URL`).
#'   \item **Assays**: `PValue`, `QValueFDRBH`, `QValueFDRBY`,
#'     `QValueBonferroni`, `GenesInTermInQuery` and `Genes`. Terms that were not
#'     enriched in a cluster are `NA`.
#'   \item **metadata**: `DEGs` (the standardized input with columns `cluster`,
#'     `gene`, `logFC`, `p_val`), `toppInput` and, after [toppFun()], `toppParams`.
#' }
#'
#' If `input_data` is the results data.frame returned by [toppFun()], it is
#' converted directly (the `*_col` arguments are ignored, and the `Cluster` and
#' optional `Direction` columns are used). DEGs can then be attached afterwards
#' with [addDEGs()].
#'
#' @param input_data A DEG data.frame, a marker list/data.frame (see the `type`
#'   argument of [toppFun()]), or a results data.frame returned by [toppFun()]
#' @param type One of "degs", "marker_df", "marker_list" or "results". By default
#'   this is detected from `input_data`.
#' @param cluster_col Column name for the groups of cells (e.g. cluster or celltype)
#' @param gene_col Column name for genes (e.g. gene or feature)
#' @param logFC_col Column name for the avg log FC column
#' @param p_val_col Column name for the p-value or adjusted p-value (preferred)
#' @param key_type Gene name format, "SYMBOL" or "ENTREZ"
#' @return A SummarizedExperiment
#' @seealso [toppFun()], [toppData()], [addDEGs()], [termGenes()]
#' @examples
#' data("ifnb.de")
#' se <- toppSE(ifnb.de,
#'     cluster_col = "celltype",
#'     gene_col = "gene",
#'     p_val_col = "p_val_adj",
#'     logFC_col = "avg_log2FC"
#' )
#' se
#'
#' # From existing toppFun() results
#' data("toppdata.ifnb")
#' se_results <- toppSE(toppdata.ifnb)
#' dim(se_results)
#' @export
toppSE <- function(input_data,
                   type = NULL,
                   cluster_col = "cluster",
                   gene_col = "gene",
                   logFC_col = "avg_log2FC",
                   p_val_col = "p_val_adj",
                   key_type = "SYMBOL") {
    .check_se_pkgs()
    if (is.null(type)) {
        type <- if (.looks_like_results(input_data)) {
            "results"
        } else if (is.data.frame(input_data)) {
            "degs"
        } else if (is.list(input_data)) {
            "marker_df"
        } else {
            "marker_list"
        }
    }
    if (!(type %in% c("degs", "marker_df", "marker_list", "results"))) {
        stop("Please ensure the parameter `type` is one of: degs, marker_df, marker_list, or results.")
    }

    md <- list(toppInput = list(type = type, key_type = key_type))

    if (type == "results") {
        if (!.looks_like_results(input_data) || !("Cluster" %in% colnames(input_data))) {
            stop("Results input must be a data.frame from toppFun() with the columns ",
                paste(c(.topp_required_cols, "Cluster"), collapse = ", "), ".",
                call. = FALSE
            )
        }
        results <- as.data.frame(input_data, stringsAsFactors = FALSE)
        direction <- if ("Direction" %in% colnames(results)) as.character(results$Direction) else rep("all", nrow(results))
        direction[is.na(direction)] <- "all"
        pairs <- unique(data.frame(
            Cluster = as.character(results$Cluster), Direction = direction,
            stringsAsFactors = FALSE
        ))
        col_df <- data.frame(
            Cluster = pairs$Cluster, Direction = pairs$Direction,
            n_input_genes = NA_integer_, stringsAsFactors = FALSE
        )
        return(.build_topp_se(results, col_df, md))
    }

    degs <- if (type == "degs") {
        .standardize_degs(input_data, cluster_col, gene_col, logFC_col, p_val_col)
    } else {
        .standardize_markers(input_data, type)
    }
    if (nrow(degs) == 0) {
        stop("No genes found in `input_data`.", call. = FALSE)
    }
    clusters <- unique(degs$cluster)
    col_df <- data.frame(
        Cluster = clusters, Direction = "all",
        n_input_genes = .count_input_genes(degs, clusters), stringsAsFactors = FALSE
    )
    md$DEGs <- degs
    .build_topp_se(NULL, col_df, md)
}

#' Extract toppFun results from a SummarizedExperiment as a data.frame
#'
#' @description
#' Converts the term-by-cluster enrichment results stored in a [toppSE()] object
#' back into the long data.frame format returned by [toppFun()]. For objects that
#' were created with [addToppData()], the data.frame stored in the metadata slot
#' is returned.
#'
#' @param se A SummarizedExperiment created with [toppSE()] or [addToppData()]
#' @param slot_name For objects created with [addToppData()], the metadata slot name
#'   containing toppData (default: "toppData")
#' @return A data.frame with one row per term and cluster
#' @seealso [toppSE()]
#' @examples
#' data("toppdata.ifnb")
#' se <- toppSE(toppdata.ifnb)
#' head(toppData(se))
#' @export
toppData <- function(se, slot_name = "toppData") {
    .check_se_pkgs()
    if (!inherits(se, "SummarizedExperiment")) {
        stop("`se` must be a SummarizedExperiment or SingleCellExperiment object")
    }
    if (.is_toppSE(se)) {
        return(.se_to_results(se))
    }

    metadata_obj <- S4Vectors::metadata(se)
    if (!slot_name %in% names(metadata_obj)) {
        stop("toppData not found in metadata slot '", slot_name, "'. ",
            "Available metadata slots: ", paste(names(metadata_obj), collapse = ", "), ". ",
            "Use addToppData() to add results or specify correct slot_name.",
            call. = FALSE
        )
    }
    extracted_data <- metadata_obj[[slot_name]]
    if (!is.data.frame(extracted_data)) {
        stop("Data in metadata slot '", slot_name, "' is not a data.frame. ",
            "Expected toppData results stored as data.frame.",
            call. = FALSE
        )
    }
    extracted_data
}

#' Add differential expression results to a toppSE object
#'
#' @description
#' Stores DEG results, with standardized column names, in a [toppSE()] object. This
#' lets the hit genes of each enriched term be linked back to their fold changes
#' and p-values with [termGenes()], and lets [toppFun()] be (re)run on the object.
#' Any DEGs already stored in the object are replaced.
#'
#' @param se A SummarizedExperiment created with [toppSE()]
#' @param degs A data.frame of differential expression results
#' @param cluster_col Column name for the groups of cells (e.g. cluster or celltype).
#'   These values must match the clusters in `se`.
#' @param gene_col Column name for genes (e.g. gene or feature)
#' @param logFC_col Column name for the avg log FC column
#' @param p_val_col Column name for the p-value or adjusted p-value (preferred)
#' @param key_type Gene name format, "SYMBOL" or "ENTREZ". `NULL` keeps the format
#'   already recorded in the object.
#' @return The SummarizedExperiment with the DEGs stored in `metadata(se)$DEGs`
#' @seealso [toppSE()], [termGenes()]
#' @examples
#' data("toppdata.ifnb")
#' data("ifnb.de")
#' se <- toppSE(toppdata.ifnb)
#' se <- addDEGs(se, ifnb.de,
#'     cluster_col = "celltype",
#'     gene_col = "gene",
#'     p_val_col = "p_val_adj",
#'     logFC_col = "avg_log2FC"
#' )
#' @export
addDEGs <- function(se,
                    degs,
                    cluster_col = "cluster",
                    gene_col = "gene",
                    logFC_col = "avg_log2FC",
                    p_val_col = "p_val_adj",
                    key_type = NULL) {
    .check_se_pkgs()
    if (!.is_toppSE(se)) {
        stop("`se` must be a SummarizedExperiment created with toppSE().", call. = FALSE)
    }
    degs <- .standardize_degs(degs, cluster_col, gene_col, logFC_col, p_val_col)

    se_clusters <- unique(as.character(SummarizedExperiment::colData(se)$Cluster))
    deg_clusters <- unique(degs$cluster)
    if (length(intersect(se_clusters, deg_clusters)) == 0) {
        stop("None of the clusters in `degs` match the clusters in `se`.\n",
            "Clusters in se: ", paste(se_clusters, collapse = ", "), "\n",
            "Clusters in degs: ", paste(deg_clusters, collapse = ", "),
            call. = FALSE
        )
    }
    missing_clusters <- setdiff(se_clusters, deg_clusters)
    if (length(missing_clusters) > 0) {
        warning("No DEGs found for cluster(s): ", paste(missing_clusters, collapse = ", "), call. = FALSE)
    }

    md <- S4Vectors::metadata(se)
    md$DEGs <- degs
    md$toppInput <- list(
        type = "degs",
        key_type = key_type %||% md$toppInput$key_type %||% "SYMBOL"
    )
    S4Vectors::metadata(se) <- md
    SummarizedExperiment::colData(se)$n_input_genes <- .count_input_genes(
        degs, as.character(SummarizedExperiment::colData(se)$Cluster)
    )
    se
}

#' Get the genes behind an enriched term, with their DEG statistics
#'
#' @description
#' Returns the genes that hit a ToppGene term in each cluster, joined with the
#' fold change and p-value from the DEGs stored in the object (see [toppSE()] and
#' [addDEGs()]). This is the starting point for gene-level visualizations of a
#' term, such as heatmaps, bar plots or volcano plots.
#'
#' @param se A SummarizedExperiment created with [toppSE()]
#' @param term The term ID or name
#' @param cluster Optional cluster name(s) to restrict the results to
#' @param category Optional category name(s), useful if a term appears in more than one
#' @return A data.frame with one row per gene, cluster and term, with the columns
#'   `Cluster`, (`Direction`), `Category`, `ID`, `Name`, `gene`, `logFC` and `p_val`.
#'   `logFC` and `p_val` are `NA` for genes without a matching DEG.
#' @details Genes are matched to the DEGs by name, and ToppGene returns gene
#'   symbols, so DEGs given as Entrez IDs will not match.
#' @seealso [toppSE()], [addDEGs()]
#' @examples
#' data("toppdata.ifnb")
#' data("ifnb.de")
#' se <- toppSE(toppdata.ifnb)
#' se <- addDEGs(se, ifnb.de,
#'     cluster_col = "celltype",
#'     gene_col = "gene",
#'     p_val_col = "p_val_adj",
#'     logFC_col = "avg_log2FC"
#' )
#' termGenes(se, "cytokine receptor binding", cluster = "CD14 Mono")
#' @export
termGenes <- function(se, term, cluster = NULL, category = NULL) {
    .check_se_pkgs()
    if (!.is_toppSE(se)) {
        stop("`se` must be a SummarizedExperiment created with toppSE().", call. = FALSE)
    }
    results <- toppData(se)
    hit <- results$ID %in% term | results$Name %in% term
    if (!is.null(category)) {
        hit <- hit & results$Category %in% category
    }
    if (!is.null(cluster)) {
        hit <- hit & results$Cluster %in% cluster
    }
    if (!any(hit)) {
        stop("No enrichment results found for the requested term/cluster.", call. = FALSE)
    }
    results <- results[hit, , drop = FALSE]

    if (!is.character(results$Genes)) {
        stop("The enrichment results have no `Genes` column; re-run toppFun() to include the genes in each term.", call. = FALSE)
    }
    gene_lists <- strsplit(results$Genes, ",\\s*")
    id_cols <- intersect(c("Cluster", "Direction", "Category", "ID", "Name"), colnames(results))
    out <- results[rep(seq_len(nrow(results)), lengths(gene_lists)), id_cols, drop = FALSE]
    out$gene <- unlist(gene_lists)
    rownames(out) <- NULL

    md <- S4Vectors::metadata(se)
    degs <- md$DEGs
    out$logFC <- NA_real_
    out$p_val <- NA_real_
    if (!is.null(degs) && nrow(degs) > 0) {
        m <- match(paste(out$Cluster, out$gene), paste(degs$cluster, degs$gene))
        out$logFC <- degs$logFC[m]
        out$p_val <- degs$p_val[m]
        if (all(is.na(m)) && identical(md$toppInput$key_type, "ENTREZ")) {
            warning("The stored DEGs use Entrez IDs but ToppGene returns gene symbols, so no genes could be matched.",
                call. = FALSE
            )
        }
    }
    out
}

#' @rdname toppFun
#' @export
toppFun.SummarizedExperiment <- function(input_data,
                                         topp_categories = NULL,
                                         direction_mode = "all",
                                         num_genes = 1000,
                                         pval_cutoff = 0.5,
                                         fc_cutoff = 0,
                                         fc_filter = "ALL",
                                         clusters = NULL,
                                         correction = "FDR",
                                         key_type = NULL,
                                         min_genes = 2,
                                         max_genes = 1500,
                                         max_results = 50,
                                         verbose = TRUE,
                                         ...) {
    .check_se_pkgs()
    se <- input_data
    if (!.is_toppSE(se)) {
        stop("`input_data` must be a SummarizedExperiment created with toppSE().", call. = FALSE)
    }
    md <- S4Vectors::metadata(se)
    degs <- md$DEGs
    if (is.null(degs) || nrow(degs) == 0) {
        stop("No DEGs or markers are stored in the object. Add them with addDEGs().", call. = FALSE)
    }
    type <- md$toppInput$type %||% "degs"
    key_type <- key_type %||% md$toppInput$key_type %||% "SYMBOL"
    if (!(direction_mode %in% c("all", "split"))) {
        stop("Please ensure the parameter `direction_mode` is one of: all, split.", call. = FALSE)
    }
    if (direction_mode == "split" && type != "degs") {
        stop("`direction_mode = \"split\"` requires DEGs with fold changes, not marker lists.", call. = FALSE)
    }

    all_clusters <- unique(degs$cluster)
    if (!is.null(clusters)) {
        clusters <- as.character(clusters)
        unknown <- setdiff(clusters, all_clusters)
        if (length(unknown) > 0) {
            stop("Cluster(s) not found in the stored DEGs: ", paste(unknown, collapse = ", "), call. = FALSE)
        }
        all_clusters <- all_clusters[all_clusters %in% clusters]
    }
    degs <- degs[degs$cluster %in% all_clusters, , drop = FALSE]

    if (nrow(se) > 0) {
        warning("Replacing the enrichment results already stored in the object.", call. = FALSE)
    }

    query_args <- list(
        topp_categories = topp_categories, direction_mode = direction_mode,
        num_genes = num_genes, pval_cutoff = pval_cutoff, fc_cutoff = fc_cutoff,
        fc_filter = fc_filter, correction = correction, key_type = key_type,
        min_genes = min_genes, max_genes = max_genes, max_results = max_results,
        verbose = verbose
    )
    results <- if (type == "degs") {
        do.call(toppFun.default, c(
            list(
                input_data = degs, type = "degs", cluster_col = "cluster", gene_col = "gene",
                p_val_col = "p_val", logFC_col = "logFC"
            ),
            query_args
        ))
    } else {
        markers <- split(degs$gene, factor(degs$cluster, levels = all_clusters))
        do.call(toppFun.default, c(list(input_data = markers, type = "marker_df"), query_args))
    }

    directions <- if (direction_mode == "split") c("up", "down") else "all"
    col_df <- data.frame(
        Cluster = rep(all_clusters, each = length(directions)),
        Direction = rep(directions, times = length(all_clusters)),
        stringsAsFactors = FALSE
    )
    col_df$n_input_genes <- .count_input_genes(degs, col_df$Cluster)

    md$toppParams <- c(
        list(timestamp = Sys.time(), type = type),
        query_args[setdiff(names(query_args), c("topp_categories", "verbose"))],
        list(topp_categories = topp_categories %||% get_ToppCats())
    )
    .build_topp_se(results, col_df, md)
}
