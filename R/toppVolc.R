## Build one volcano plot for a single gene set
.volcano_one <- function(degs, genes, title, subtitle, logFC_threshold,
                         pval_threshold, label, label_max, point_alpha) {
    degs <- degs[!is.na(degs$logFC) & !is.na(degs$p_val), , drop = FALSE]
    degs$Highlight <- ifelse(degs$gene %in% genes, "In gene set", "Other")
    degs$Highlight <- factor(degs$Highlight, levels = c("Other", "In gene set"))
    degs$neglog_p <- -log10(pmax(degs$p_val, .Machine$double.xmin))
    # draw highlighted genes on top
    degs <- degs[order(degs$Highlight), , drop = FALSE]
    highlighted <- degs[degs$Highlight == "In gene set", , drop = FALSE]

    p <- ggplot2::ggplot(degs, ggplot2::aes(x = logFC, y = neglog_p, color = Highlight)) +
        ggplot2::geom_point(alpha = point_alpha) +
        ggplot2::scale_color_manual(
            values = c("Other" = "grey70", "In gene set" = "red"),
            drop = FALSE
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(pval_threshold), linetype = "dashed", color = "grey40"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-logFC_threshold, logFC_threshold),
            linetype = "dashed", color = "grey40"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = title, subtitle = subtitle,
            x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value"
        ) +
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_text(size = 10)
        )

    if (label && nrow(highlighted) > 0) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_label_repel(
                data = highlighted,
                ggplot2::aes(label = gene),
                size = 3, max.overlaps = label_max, show.legend = FALSE
            )
        } else {
            warning("Install the 'ggrepel' package to label the highlighted genes.", call. = FALSE)
        }
    }
    p
}

#' Volcano plots highlighting the genes of ToppGene terms
#'
#' @description
#' Draws a volcano plot of the differentially expressed genes of one cluster and
#' highlights the genes that belong to an enriched ToppGene term. By default
#' one plot is made for each of the top terms in a category.
#'
#' @param toppData A SummarizedExperiment created with [toppSE()] (the DEGs
#'   are read from the object), or a data.frame of results from [toppFun()]
#'   (in which case `degs` is required). The results must include the `Genes`
#'   column.
#' @param degs A data.frame of differentially expressed genes. Only needed when
#'   `toppData` is a data.frame, or to override the DEGs stored in a toppSE.
#' @param category ToppGene category of the terms to plot (see [get_ToppCats()])
#' @param cluster The single cluster to plot, as named in the results
#' @param term Optional term name(s) or ID(s). If `NULL`, the `num_terms` most
#'   significant terms of `category` in `cluster` are used.
#' @param num_terms Number of top terms to plot when `term` is `NULL`
#' @param p_val_adj Results column used to rank the terms (default `"QValueFDRBH"`)
#' @param deg_cluster_col,gene_col,logFC_col,p_val_col Column names in `degs`.
#'   Ignored for DEGs stored in a toppSE. If `deg_cluster_col` is `NULL`, `degs`
#'   is assumed to contain a single cluster.
#' @param logFC_threshold,pval_threshold Thresholds drawn as dashed lines
#' @param label Logical, label the highlighted genes (needs the 'ggrepel' package)
#' @param label_max Maximum number of label overlaps passed to ggrepel
#' @param point_alpha Transparency of the points
#' @param combine Logical, return multiple plots as one patchwork object
#' @param ncols Number of columns when `combine = TRUE`
#' @param save Logical, save the plot(s)
#' @param save_dir Directory to save into
#' @param file_prefix Prefix of the saved file names (default `"toppVolc"`)
#' @param width,height Size of each saved plot in inches
#' @param slot_name Metadata slot holding results, for objects not created with
#'   [toppSE()]
#' @param ... Not used
#' @return A ggplot object when a single term is plotted, otherwise a named
#'   list of ggplot objects, or a patchwork object if `combine = TRUE`.
#' @examples
#' data("toppdata.ifnb")
#' data("ifnb.de")
#' toppVolc(
#'     toppdata.ifnb,
#'     degs = ifnb.de,
#'     category = "GeneOntologyMolecularFunction",
#'     cluster = "NK",
#'     num_terms = 1,
#'     deg_cluster_col = "celltype",
#'     logFC_col = "avg_log2FC",
#'     p_val_col = "p_val_adj"
#' )
#' @seealso [termGenes()], [toppPlot()], [toppSE()]
#' @export
toppVolc <- function(toppData,
                     degs = NULL,
                     category = NULL,
                     cluster = NULL,
                     term = NULL,
                     num_terms = 5,
                     p_val_adj = "QValueFDRBH",
                     deg_cluster_col = "cluster",
                     gene_col = "gene",
                     logFC_col = "avg_log2FC",
                     p_val_col = "p_val_adj",
                     logFC_threshold = 0.25,
                     pval_threshold = 0.05,
                     label = TRUE,
                     label_max = 10,
                     point_alpha = 0.6,
                     combine = FALSE,
                     ncols = 2,
                     save = FALSE,
                     save_dir = tempdir(),
                     file_prefix = NULL,
                     width = 6,
                     height = 5,
                     slot_name = "toppData",
                     ...) {
    if (is.null(category) && is.null(term)) {
        stop("Please provide a `category` (or specific `term`s) to plot.", call. = FALSE)
    }
    if (is.null(cluster) || length(cluster) != 1) {
        stop("Please provide exactly one `cluster` to plot.", call. = FALSE)
    }

    # standardized DEGs: cluster, gene, logFC, p_val
    if (is.null(degs) && inherits(toppData, "SummarizedExperiment") && .is_toppSE(toppData)) {
        degs_std <- S4Vectors::metadata(toppData)$DEGs
        if (is.null(degs_std) || nrow(degs_std) == 0) {
            stop("No DEGs are stored in the object. Add them with addDEGs() or supply `degs`.", call. = FALSE)
        }
    } else if (is.null(degs)) {
        stop("`degs` is required when `toppData` is not a toppSE with stored DEGs.", call. = FALSE)
    } else {
        degs <- as.data.frame(degs)
        if (is.null(deg_cluster_col)) {
            degs$.cluster <- cluster
            deg_cluster_col <- ".cluster"
        }
        degs_std <- .standardize_degs(degs, deg_cluster_col, gene_col, logFC_col, p_val_col)
    }

    results <- if (inherits(toppData, "SummarizedExperiment")) {
        toppData(toppData, slot_name = slot_name)
    } else if (is.data.frame(toppData)) {
        toppData
    } else {
        stop("`toppData` must be a data.frame or a SummarizedExperiment.", call. = FALSE)
    }
    if (!is.character(results$Genes)) {
        stop("The enrichment results have no `Genes` column; re-run toppFun() to include the genes in each term.", call. = FALSE)
    }
    if (!is.null(category)) {
        bad <- setdiff(category, results$Category)
        if (length(bad) > 0) {
            stop("Category not found in the results: ", paste(bad, collapse = ", "), call. = FALSE)
        }
        results <- results[results$Category %in% category, , drop = FALSE]
    }
    if (!cluster %in% results$Cluster) {
        stop("Cluster '", cluster, "' not found in the results.", call. = FALSE)
    }
    results <- results[results$Cluster == cluster, , drop = FALSE]

    if (is.null(term)) {
        results <- results[order(results[[p_val_adj]]), , drop = FALSE]
        results <- utils::head(results, num_terms)
    } else {
        results <- results[results$Name %in% term | results$ID %in% term, , drop = FALSE]
        if (nrow(results) == 0) {
            stop("None of the requested terms were found for cluster '", cluster, "'.", call. = FALSE)
        }
    }

    cluster_degs <- degs_std[degs_std$cluster == cluster, , drop = FALSE]
    if (nrow(cluster_degs) == 0) {
        stop("No DEGs found for cluster '", cluster, "'.", call. = FALSE)
    }

    plots <- lapply(seq_len(nrow(results)), function(i) {
        genes <- strsplit(results$Genes[i], ",\\s*")[[1]]
        .volcano_one(
            cluster_degs, genes,
            title = paste("ToppVolcano:", cluster),
            subtitle = paste(results$Category[i], results$Name[i], sep = " - "),
            logFC_threshold = logFC_threshold, pval_threshold = pval_threshold,
            label = label, label_max = label_max, point_alpha = point_alpha
        )
    })
    names(plots) <- results$Name

    if (isTRUE(save)) {
        prefix <- if (is.null(file_prefix)) "toppVolc" else file_prefix
        for (i in seq_along(plots)) {
            safe <- gsub("[^A-Za-z0-9]+", "_", results$Name[i])
            ggplot2::ggsave(
                filename = file.path(save_dir, paste0(prefix, "_", cluster, "_", safe, ".pdf")),
                plot = plots[[i]], width = width, height = height
            )
        }
    }

    if (length(plots) == 1) {
        return(plots[[1]])
    }
    if (isTRUE(combine)) {
        return(patchwork::wrap_plots(plots, ncol = ncols))
    }
    plots
}
