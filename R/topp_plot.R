#' Create a dotplot from toppdata results
#' NEW: this accepts only 1 category
#'
#' @param toppData A toppData results dataframe
#' @param category The topp categories to plot
#' @param clusters The cluster(s) to plot
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
#' toppPlot(toppData,
#'     category="GeneOntologyMolecularFunction",
#'     clusters=0,
#'     save=TRUE,
#'     file_prefix="MF_cluster0")
#'
#' @export
toppPlot <- function (toppData,
                      category,
                      clusters = NULL,
                      num_terms = 10,
                      p_val_adj="BH",
                      p_val_display="log",
                      save = FALSE,
                      save_dir = NULL,
                      width = 5,
                      height = 6,
                      file_prefix = NULL,
                      y_axis_text_size=8,
                      combine = FALSE,
                      ncols = NULL) {
  #check for correct columns
  test_cols = c("Category", "Name", "PValue", "GenesInTerm", "GenesInQuery", "GenesInTermInQuery")
  for (t in test_cols) {
    if (!(t %in% colnames(toppData))) {
      stop(paste0("The column ", t, " is missing from the toppData, please correct and retry."))
    }
  }
  #parse category - make sure it's in data and there is only 1
  if (!(category %in% unique(toppData$Category))){
    stop(paste0("Category ", category, " not found in the data. Please select one of ",
                stringr::str_c(unique(toppData$Category), collapse = ", ")))
  } else if (is.null(category)){
    stop(paste0("Please select one of these categories: ",
                stringr::str_c(unique(toppData$Category), collapse = ", ")))
  }

  #parse clusters
  if (is.null(clusters)) {
    if ("Cluster" %in% colnames(toppData)) {
      clusters = unique(toppData$Cluster)
    }
  }

  #parse pvalue

  tmp_data <- toppData

  if (!(p_val_adj %in% c("BH", "Bonferroni", "BY", "none", "None", "log"))) {
    cat("P value adjustment not found - using 'BH' by default. For no adjustment, use p_val_adj = 'none'.")
  }


  p_val_col = switch(p_val_adj,
                     "BH" = "QValueFDRBH",
                     "Bonferroni" = "QValueBonferroni",
                     "BY" = "QvalueFDRBY",
                     "none" = "PValue",
                     "None" = "PValue",
                     "QValueFDRBH")


  if (p_val_display == 'log') {
    color_label = "-Log10(FDR)"
    tmp_data <- tmp_data |>
      dplyr::mutate(nlog10_fdr = -log10(!!as.name(p_val_col)))
    p_val_display_column = "nlog10_fdr"
  } else if (p_val_col %in% c("QValueFDRBH", "QValueBonferroni", "QvalueFDRBY")) {
    color_label = "Adj. P-value"
    p_val_display_column = p_val_col
  } else {
    color_label = "P-value"
    p_val_display_column = p_val_col
  }

  #parse save
  if (isTRUE(save)) {
    if (is.null(save_dir)) {
      output_dir = getwd()
    } else {
      output_dir = save_dir
    }
  }

  ###
  ###MAIN PLOTTING SECTION
  ###
  if (length(clusters) > 1) {
    if (!isTRUE(combine)) {
      cat("Multiple clusters entered: function returns a list of ggplots\n")
    }
    overall_plot_list = list()
      for (c in clusters) {
        overall_plot_list[[c]] <- tmp_data |>
          dplyr::filter(Cluster == c) |>
          dplyr::filter(Category == category) |>
          dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
          dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
          dplyr::arrange(-geneRatio) |>
          utils::head(num_terms) |>

          ggplot2::ggplot(mapping = aes(
            x = geneRatio,
            y = Name
          )) +
          ggplot2::geom_segment(aes(xend=0, yend=Name)) +
          ggplot2::geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_display_column))) +
          viridis::scale_color_viridis(option = "C")  +
          ggplot2::theme_bw() +
          ggplot2::ylab(category) +
          ggplot2::ggtitle(stringr::str_glue("Cluster {c}")) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_axis_text_size)) +
          ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = FALSE)) +
          ggplot2::labs(color=color_label, size = "Genes from Query\n in Gene Set")

        if (isTRUE(save)) {
          if (is.null(file_prefix)) {
            save_filename = stringr::str_glue("{category}_{c}_toppDotPlot.pdf")
          } else {
            save_filename = stringr::str_glue("{file_prefix}_{c}.pdf")
          }
          ggplot2::ggsave(filename = file.path(output_dir, save_filename),
                          width = width, height=height)
        }
      }

    if (isTRUE(combine)){
      if (is.null(ncols)) {
        ncols = min(3, length(overall_plot_list))
      }
      combined_plots = patchwork::wrap_plots(overall_plot_list, ncol = ncols) +
        patchwork::plot_annotation(title=category)
      return (combined_plots)
    } else {
      return (overall_plot_list)
    }
  } else if (length(clusters) == 1) {
    c = clusters[1]
      single_plot <- tmp_data |>
        dplyr::filter(Cluster == c) |>
        dplyr::filter(Category == category) |>
        dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
        dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
        dplyr::arrange(-geneRatio) |>
        utils::head(num_terms) |>

        ggplot2::ggplot(mapping = aes(
          x = geneRatio,
          y = Name
        )) +
        ggplot2::geom_segment(aes(xend=0, yend=Name)) +
        ggplot2::geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_display_column))) +
        viridis::scale_color_viridis(option = "C")  +
        ggplot2::theme_bw() +
        ggplot2::ylab(category) +
        ggplot2::ggtitle(stringr::str_glue("Cluster: {c}")) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = y_axis_text_size)) +
        ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = FALSE)) +
        ggplot2::labs(color=color_label, size = "Genes from Query\n in Gene Set")

      if (isTRUE(save)) {
        if (is.null(file_prefix)) {
          save_filename = stringr::str_glue("{category}_{c}_toppDotPlot.pdf")
        } else {
          save_filename = stringr::str_glue("{file_prefix}_{category}_{c}_toppDotPlot.pdf")
        }
        ggplot2::ggsave(filename = file.path(output_dir, save_filename),
                        width = width, height=height)
        }
    return (single_plot)
    # if (length(category_plot_list) == 1){
    #   category_plot_list[[1]]
    # } else {
    #   return (category_plot_list)
    # }
  }
}

#' Create a balloon plot from toppdata results
#'
#' @param toppData A toppData results dataframe
#' @param categories The topp categories to plot
#' @param balloons Number of balloons per group to plot
#' @param filename Filename of the saved balloon plot
#' @param save Save the balloon plot if TRUE
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
#' toppBalloon(toppData, balloons = 3, save = TRUE, filename = "Balloon_plot")
#'
#' @export
toppBalloon <- function (toppData,
                         categories = NULL,
                         balloons = 3,
                         x_axis_text_size = 6,
                         filename = NULL,
                         save = FALSE,
                         height = 5,
                         width = 10) {

  if (is.null(categories)) {
    categories <- unique(toppData[["Category"]])
  }
  balloon_list <- list()
  for (cat in categories) {
    cat("Balloon Plot:", cat)
    balloon_list[[cat]] <- toppData |>
      dplyr::filter(Category == cat) |>
      dplyr::mutate(nlog10_fdr = -log10(QValueFDRBH)) |>
      dplyr::group_by(Cluster) |>
      dplyr::slice_max(order_by=nlog10_fdr, n=balloons, with_ties = FALSE) |>
      dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
      ggplot(aes(
        x=forcats::fct_reorder(Name, Cluster),
        y=Cluster,
      )) +
      ggplot2::geom_point(aes(size=geneRatio, color=nlog10_fdr)) +
      viridis::scale_color_viridis(option = "C") +
      ggplot2::labs(color="-Log10(FDR)", size="Gene Ratio") +
      ggplot2::xlab(cat) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30, whitespace_only = FALSE)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = x_axis_text_size, angle = 60, hjust=1.1, vjust=1.05),
                     panel.border = ggplot2::element_rect(color = NA))

      if (is.null(filename)) {
        plot_filename <- stringr::str_glue("toppBalloon_{cat}.pdf")
      } else {
        plot_filename = stringr::str_glue("{filename}_{cat}.pdf")
      }
      if (isTRUE(save)){
        ggplot2::ggsave(plot_filename, height = height, width = width)
        cat(" saved to:", plot_filename)
      }
    cat("\n")

  }
  # if (length(balloon) == 1) {
  #   print(balloon[1])
  # }
  #return (balloon)
  if (length(balloon_list) == 1){
    balloon_list[[1]]
  } else {
    return (balloon_list)
  }
}
