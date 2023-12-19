#' Create a dotplot from toppdata results
#'
#' @param toppData A toppData results dataframe
#' @param category The topp category to plot
#' @param clusters The cluster(s) to plot
#' @param p_val_adj The P-value correction method: "BH", "Bonferroni", "BY", or "none"
#' @param num_terms The number of terms from the toppData results to be plotted, per cluster
#' @return ggplot
#' @importFrom ggplot2 ggplot geom_point geom_segment scale_color_gradient theme_bw scale_y_discrete labs aes
#' @importFrom dplyr filter mutate
#' @importFrom stringr str_wrap
#' @importFrom forcats fct_reorder
#' @importFrom viridis scale_color_viridis
#' @export
toppPlot <- function (toppData,
                      category = NULL,
                      clusters = NULL,
                      num_terms = 10,
                      p_val_adj="BH",
                      save = FALSE,
                      save_dir = NULL,
                      width = 7) {
  #check toppData object - check for columns? - check that clusters and category are in data
  test_cols = c("Category", "Name", "PValue", "GenesInTerm", "GenesInQuery", "GenesInTermInQuery")
  for (t in test_cols) {
    if (!(t %in% colnames(toppData))) {
      stop(paste0("The column ", t, " is missing from the toppData, please correct and retry."))
    }
  }
  #parse category
  if (!(is.null(category))){
    if (length(c(category)) > 1) {
      stop("Please select one ToppFun category to plot (e.g. category = 'GeneOntologyMolecularFunction').")
    }
    else if (length(c(category)) == 1) {
      category = category[1]
    }
    else {
      stop("Invalid input data - please check your toppData dataframe")
    }
  }

  #parse pvalue

  if (!(p_val_adj %in% c("BH", "Bonferroni", "BY", "none", "None", "log"))) {
    cat("P value adjustment not found - using 'BH' by default. For no adjustment, use p_val_adj = 'none'.")
  }

  if ("nlog10_fdr" %in% colnames(toppData)) {
    p_val_adj = "log"
  }

  p_val_col = switch(p_val_adj,
                     "BH" = "QValueFDRBH",
                     "Bonferroni" = "QValueBonferroni",
                     "BY" = "QvalueFDRBY",
                     "none" = "PValue",
                     "None" = "PValue",
                     "log" = "nlog10_fdr",
                     "QValueFDRBH")

  if (p_val_col %in% c("QValueFDRBH", "QValueBonferroni", "QvalueFDRBY")) {
    color_label = "Adj. P-value"
  } else if (p_val_col == "nlog10_fdr") {
    color_label = "-Log10(FDR)"
  } else {
    color_label = "P-value"
  }

  #parse save
  if (isTRUE(save)) {
    if (is.null(save_dir)) {
      output_dir = getwd()
    } else {
      output_dir = save_dir
    }
  }

  #parse clusters
  if (is.null(clusters)) {
    if ("Cluster" %in% colnames(toppData)) {
      clusters = unique(toppData$Cluster)
    }
  }

  clusters <- clusters |> as.character()

  if (length(clusters) > 1) {
    cat("Multiple clusters entered: function returns a list of ggplots")

    plot_list = list()
    for (cat in category) {
      for (c in clusters) {
        plot_list[[c]] <- toppData |>
          dplyr::filter(Cluster == c) |>
          dplyr::filter(Category == category) |>
          dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
          dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
          dplyr::arrange(-geneRatio) |>
          head(num_terms) |>

          ggplot2::ggplot(mapping = aes(
            x = geneRatio,
            y = Name
          )) +
          ggplot2::geom_segment(aes(xend=0, yend=Name)) +
          ggplot2::geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_col))) +
          viridis::scale_color_viridis()  +
          ggplot2::theme_bw() +
          ggplot2::ylab(cat) +
          ggplot2::ggtitle(stringr::str_glue("Cluster {c}")) +
          ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4)) +
          ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = F)) +
          ggplot2::labs(color=color_label, size = "Genes from Query\n in Gene Set")

        if (isTRUE(save)) {
          filename = stringr::str_glue("{cat}_{c}_toppDotPlot.pdf")
          ggplot2::ggsave(filename = file.path(output_dir, filename),
                          width = width)
        }
      }
    }
    return (plot_list)
  } else {
      for (cat in category) {
        for (c in clusters) {
          dotplot <- toppData |>
            dplyr::filter(Cluster == c) |>
            dplyr::filter(Category == category) |>
            dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
            dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
            dplyr::arrange(-geneRatio) |>
            head(num_terms) |>

            ggplot2::ggplot(mapping = aes(
              x = geneRatio,
              y = Name
            )) +
            ggplot2::geom_segment(aes(xend=0, yend=Name)) +
            ggplot2::geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_col))) +
            viridis::scale_color_viridis()  +
            ggplot2::theme_bw() +
            ggplot2::ylab(cat) +
            ggplot2::ggtitle(stringr::str_glue("Cluster {c}")) +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4)) +
            ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = F)) +
            ggplot2::labs(color=color_label, size = "Genes from Query\n in Gene Set")

          if (isTRUE(save)) {
            filename = stringr::str_glue("{cat}_{c}_toppDotPlot.pdf")
            ggplot2::ggsave(filename = file.path(output_dir, filename),
                            width = width)
          }
        }
      }
    return (dotplot)
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
#' @return ggplot
#' @importFrom ggplot2 ggplot geom_point theme scale_x_discrete theme_bw labs aes xlab ggsave
#' @importFrom dplyr group_by mutate slice_max
#' @importFrom stringr str_wrap str_glue
#' @importFrom forcats fct_reorder
#' @importFrom viridis scale_color_viridis
#' @export
toppBalloon <- function (toppData,
                         categories = NULL,
                         balloons = 3,
                         filename = NULL,
                         save = TRUE,
                         height = 5,
                         width = 10) {

  if (is.null(categories)) {
    categories <- unique(toppData[["Category"]])
  }
  balloon <- list()
  for (cat in categories) {
    cat("Balloon Plot:", cat)
  balloon[[cat]] <- toppData |>
    dplyr::filter(Category == cat) |>
    dplyr::group_by(Cluster) |>
    dplyr::slice_max(order_by=-nlog10_fdr, n=balloons, with_ties = F) |>
    dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
    ggplot(aes(
      x=forcats::fct_reorder(Name, Cluster),
      y=Cluster,
    )) +
    ggplot2::geom_point(aes(size=geneRatio, color=nlog10_fdr)) +
    viridis::scale_color_viridis(option = "C") +
    ggplot2::labs(color="FDR", size="Genes") +
    ggplot2::xlab(cat) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30, whitespace_only = F)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 6, angle = 60, hjust=1.1, vjust=1.05),
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
  if (length(balloon) == 1) {
    print(balloon[1])
  }
  return (balloon)
}


# else if ("Cluster" %in% colnames(toppData)) {
#     toppData |>
#       dplyr::filter(Cluster == clusters) |>
#       dplyr::filter(Category == category) |>
#       dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
#       dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
#       ggplot(mapping = aes(
#         x = geneRatio,
#         y = Name
#       )) +
#       geom_segment(aes(xend=0, yend=Name)) +
#       geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_col))) +
#       viridis::scale_color_viridis()  +
#       theme_bw() +
#       scale_y_discrete(labels = function(x) str_wrap(x, width = 20, whitespace_only = F)) +
#       labs(color=color_label, size = "Genes from Query\n in Gene Set")
#     if (isTRUE(save)) {
#       filename = str_glue("{cat}_{c}_toppDotPlot.pdf")
#       ggplot2::ggsave(filename = file.path(output_dir, filename),
#                       width = width)
#     }
#   } else {
#     toppData |>
#       dplyr::filter(Category == category) |>
#       dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
#       dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
#       ggplot(mapping = aes(
#         x = geneRatio,
#         y = Name
#       )) +
#       geom_segment(aes(xend=0, yend=Name)) +
#       geom_point(mapping = aes(size=GenesInTermInQuery, color=!!as.name(p_val_col))) +
#       viridis::scale_color_viridis() +
#       theme_bw() +
#       scale_y_discrete(labels = function(x) str_wrap(x, width = 20, whitespace_only = F)) +
#       labs(color=color_label, size = "Genes from Query in Gene Set")
#     if (isTRUE(save)) {
#       filename = str_glue("{category}_{unique}_toppDotPlot.pdf")
#       ggplot2::ggsave(filename = file.path(output_dir, filename),
#                       width = width)
#     }
#   }
#
# }
#
#
