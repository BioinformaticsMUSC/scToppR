## helper function to create dotplot for a single cluster
.sctoppr_dotplot <- function(
  topp_results,
  category,
  groupby_col,
  cluster = NULL,
  num_terms = 10,
  p_val_display_column = "FDR_BH",
  color_label = "Adjusted P-value",
  y_axis_text_size = 10
) {
  return_plot <- topp_results |>
    dplyr::filter(!!as.name(groupby_col) == cluster) |>
    dplyr::filter(Category == category) |>
    dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
    dplyr::arrange(-!!as.name(p_val_display_column)) |>
    utils::head(num_terms) |>
    ggplot2::ggplot(mapping = aes(
      x = geneRatio,
      y = forcats::fct_reorder(Name, geneRatio)
    )) +
    ggplot2::geom_segment(aes(xend = 0, yend = Name)) +
    ggplot2::geom_point(mapping = aes(size = GenesInTermInQuery, color = !!as.name(p_val_display_column))) +
    viridis::scale_color_viridis(option = "C") +
    ggplot2::theme_bw() +
    ggplot2::ylab(category) +
    ggplot2::ggtitle(stringr::str_glue("Cluster {cluster}")) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_axis_text_size)) +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = FALSE)) +
    ggplot2::labs(color = color_label, size = "Genes from Query\n in Gene Set")

  return(return_plot)
}

# helper function to save plot
.sctoppr_saveplot <- function(
  plot,
  category,
  save_dir,
  file_prefix = "toppPlot",
  width = 8,
  height = 6,
  verbose = TRUE
) {
  plot_filename <- stringr::str_glue("{file_prefix}_{category}.pdf")
  message("plot_filename:", plot_filename)
  print(save_dir)
  print(file.path(save_dir, plot_filename))
  ggplot2::ggsave(
    filename = file.path(save_dir, plot_filename),
    plot = plot,
    width = width, height = height
  )
  if (isTRUE(verbose)) {
    save_msg <- stringr::str_glue("Plot saved to {file.path(save_dir, plot_filename)}\n")
    cat(save_msg)
  }
}


#' Create a dotplot from toppdata results
#' NEW: this accepts only 1 category
#'
#' @param toppData A toppData results dataframe
#' @param category The topp categories to plot
#' @param clusters The cluster(s) to plot
#' @param cluster_col The column name for clusters (default: "Cluster")
#' @param p_val_adj The P-value correction method: "BH", "Bonferroni", "BY", or "none"
#' @param p_val_display If "log", display the p-value in terms of -log10(p_value)
#' @param num_terms The number of terms from the toppData results to be plotted, per cluster
#' @param save Whether to save the file automatically
#' @param save_dir Directory to save file
#' @param width width of the saved file (inches)
#' @param height height of the saved file (inches)
#' @param file_prefix file prefix if saving the plot - the cluster name is also added automatically
#' @param combine If TRUE and multiple clusters selected, return a patchwork object of all plots; if FALSE return list of plots
#' @param ncols If patchwork element returned, number of columns for subplots
#' @param y_axis_text_size Size of the Y axis text - for certain categories, it's helpful to decrease this
#' @return ggplot
#' @importFrom ggplot2 ggplot geom_point geom_segment scale_color_gradient theme_bw scale_y_discrete labs aes
#' @importFrom dplyr filter mutate
#' @importFrom stringr str_wrap
#' @importFrom forcats fct_reorder
#' @importFrom viridis scale_color_viridis
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom utils head
#' @examples
#' data("toppData")
#' \dontrun{
#' toppPlot(toppData,
#'   category = "GeneOntologyMolecularFunction",
#'   clusters = 0,
#'   save = TRUE,
#'   file_prefix = "MF_cluster0",
#'   save_dir = "/path/to/save/dir"
#' )
#' }
#' @export
toppPlot <- function(toppData,
                     category,
                     clusters = NULL,
                     cluster_col = "Cluster",
                     num_terms = 10,
                     p_val_adj = "BH",
                     p_val_display = "log",
                     save = FALSE,
                     save_dir = NULL,
                     width = 5,
                     height = 6,
                     file_prefix = "toppPlot",
                     y_axis_text_size = 8,
                     combine = FALSE,
                     ncols = NULL) {
  # check for correct columns
  test_cols <- c("Category", "Name", "PValue", "GenesInTerm", "GenesInQuery", "GenesInTermInQuery")
  for (t in test_cols) {
    if (!(t %in% colnames(toppData))) {
      stop(paste0("The column ", t, " is missing from the toppData, please correct and retry."))
    }
  }
  GROUPBY_COL <- cluster_col
  if (!(GROUPBY_COL %in% colnames(toppData))) {
    stop("Invalid cluster column: ", GROUPBY_COL, " not found in toppData. Select an existing column with `cluster_col = `")
  }
  # parse category - make sure it's in data and there is only 1
  if (!(category %in% unique(toppData$Category))) {
    stop(paste0(
      "Category ", category, " not found in the data. Please select one of ",
      stringr::str_c(unique(toppData$Category), collapse = ", ")
    ))
  } else if (is.null(category)) {
    stop(paste0(
      "Please select one of these categories: ",
      stringr::str_c(unique(toppData$Category), collapse = ", ")
    ))
  }

  # parse clusters - if null, then get all unique clusters
  if (is.null(clusters)) {
    clusters <- unique(toppData[[GROUPBY_COL]])
  }
  # parse save
  if (isTRUE(save)) {
    if (is.null(save_dir)) {
      stop("Please provide a directory to save the plot using the `save_dir` parameter.")
    }
  }


  # parse pvalue
  if (!(p_val_adj %in% c("BH", "Bonferroni", "BY", "none", "None", "log"))) {
    cat("P value adjustment not found - using 'BH' by default. For no adjustment, use p_val_adj = 'none'.")
  }
  p_val_col <- switch(p_val_adj,
    "BH" = "QValueFDRBH",
    "Bonferroni" = "QValueBonferroni",
    "BY" = "QvalueFDRBY",
    "none" = "PValue",
    "None" = "PValue",
    "QValueFDRBH"
  )


  plot_data <- toppData

  if (p_val_display == "log") {
    color_label <- "-Log10(FDR)"
    plot_data <- plot_data |>
      dplyr::mutate(nlog10_fdr = -log10(!!as.name(p_val_col)))
    p_val_display_column <- "nlog10_fdr"
  } else if (p_val_col %in% c("QValueFDRBH", "QValueBonferroni", "QvalueFDRBY")) {
    color_label <- "Adj. P-value"
    p_val_display_column <- p_val_col
  } else {
    color_label <- "P-value"
    p_val_display_column <- p_val_col
  }

  ###
  ### MAIN PLOTTING SECTION
  ###

  if (length(clusters) > 1) {
    if (!isTRUE(combine)) {
      cat("Multiple clusters entered: function returns a list of ggplots\n")
    }
    overall_plot_list <- list()
    for (clust in clusters) {
      cur_plot <- .sctoppr_dotplot(
        topp_results = plot_data,
        category = category,
        groupby_col = GROUPBY_COL,
        cluster = clust,
        num_terms = num_terms,
        p_val_display_column = p_val_display_column,
        color_label = color_label,
        y_axis_text_size = y_axis_text_size
      )
      # Save individual plots if specified
      if (isTRUE(save)) {
        # add cluster name to filename prefix
        filename_with_clust <- stringr::str_glue("{file_prefix %||% 'toppPlot'}_{clust}")
        message("Saving plot for cluster:", clust)
        message("filename_with_clust:", filename_with_clust)
        .sctoppr_saveplot(
          plot = cur_plot,
          category = category,
          file_prefix = filename_with_clust,
          save_dir = save_dir,
          width = width,
          height = height
        )
      }
      overall_plot_list[[clust]] <- cur_plot
    }
    # Combine plots if specified
    if (isTRUE(combine)) {
      if (is.null(ncols)) {
        ncols <- min(3, length(overall_plot_list))
      }
      combined_plots <- patchwork::wrap_plots(overall_plot_list, ncol = ncols) +
        patchwork::plot_annotation(title = category)
      return(combined_plots)
    } else {
      return(overall_plot_list)
    }
    # only one cluster
  } else if (length(clusters) == 1) {
    c <- clusters[1]
    single_plot <- .sctoppr_dotplot(
      topp_results = plot_data,
      category = category,
      groupby_col = GROUPBY_COL,
      cluster = c,
      num_terms = num_terms,
      p_val_display_column = p_val_display_column,
      color_label = color_label,
      y_axis_text_size = y_axis_text_size
    )

    # Save individual plot if specified
    if (isTRUE(save)) {
      .sctoppr_saveplot(
        plot = single_plot,
        category = category,
        file_prefix = file_prefix,
        save_dir = save_dir,
        width = width,
        height = height
      )
    }
    return(single_plot)
  }
}

#' Create a balloon plot from toppdata results
#'
#' @param toppData A toppData results dataframe
#' @param categories The topp categories to plot
#' @param cluster_col The column name for clusters (default: "Cluster")
#' @param balloons Number of balloons per group to plot
#' @param filename Filename of the saved balloon plot
#' @param save Save the balloon plot if TRUE
#' @param save_dir Directory to save the balloon plot
#' @param height Height of the saved balloon plot
#' @param width Width of the saved balloon plot
#' @param x_axis_text_size Size of the text on the x axis
#' @return ggplot
#' @importFrom ggplot2 ggplot geom_point theme scale_x_discrete theme_bw labs aes xlab ggsave
#' @importFrom dplyr group_by mutate slice_max
#' @importFrom stringr str_wrap str_glue
#' @importFrom forcats fct_reorder
#' @importFrom viridis scale_color_viridis
#' @examples
#' data("toppData")
#' \dontrun{
#' toppBalloon(toppData, balloons = 3, save = TRUE, filename = "Balloon_plot")
#' }
#'
#' @export
toppBalloon <- function(toppData,
                        categories = NULL,
                        balloons = 3,
                        x_axis_text_size = 6,
                        cluster_col = "Cluster",
                        filename = "toppBalloon",
                        save = FALSE,
                        save_dir = NULL,
                        height = 5,
                        width = 10) {
  if (is.null(categories)) {
    categories <- unique(toppData[["Category"]])
  }
  GROUPBY_COL <- cluster_col
  if (!(GROUPBY_COL %in% colnames(toppData))) {
    stop("Invalid cluster column: ", GROUPBY_COL, " not found in toppData. Select an existing column with `cluster_col = `")
  }
  # parse save
  if (isTRUE(save)) {
    if (is.null(save_dir)) {
      stop("Please provide a directory to save the plot using the `save_dir` parameter.")
    }
  }
  balloon_list <- list()
  for (cat in categories) {
    message("Creating Balloon Plot:", cat)
    balloon_list[[cat]] <- toppData |>
      dplyr::filter(Category == cat) |>
      dplyr::mutate(nlog10_fdr = -log10(QValueFDRBH)) |>
      dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
      dplyr::group_by(!!as.name(GROUPBY_COL)) |>
      dplyr::slice_max(order_by = nlog10_fdr, n = balloons, with_ties = FALSE) |>
      dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
      ggplot(aes(
        x = forcats::fct_reorder(Name, as.numeric(as.factor(!!as.name(GROUPBY_COL)))),
        y = !!as.name(GROUPBY_COL),
      )) +
      ggplot2::geom_point(aes(size = geneRatio, color = nlog10_fdr)) +
      viridis::scale_color_viridis(option = "C") +
      ggplot2::labs(color = "-Log10(FDR)", size = "Gene Ratio") +
      ggplot2::xlab(cat) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30, whitespace_only = FALSE)) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = x_axis_text_size, angle = 60, hjust = 1.1, vjust = 1.05),
        panel.border = ggplot2::element_rect(color = NA)
      )

    # Save balloon plot if specified
    if (isTRUE(save)) {
      plot_filename <- stringr::str_glue("{filename}_{cat}.pdf")
      ggplot2::ggsave(plot_filename, height = height, width = width)
      message("File saved to:", plot_filename)
    }
    cat("\n")
  }
  if (length(balloon_list) == 1) {
    balloon_list[[1]]
  } else {
    return(balloon_list)
  }
}
