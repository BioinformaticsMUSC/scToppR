#' Create a dotplot from toppdata results
#'
#' @param toppData A toppData results dataframe
#' @param category The topp category to plot
#' @param clusters The cluster(s) to plot
#' @param p_val_adj The P-value correction method: "BH", "Bonferroni", "BY", or "none"
#' @return ggplot
#' @importFrom ggplot2 ggplot geom_point geom_segment scale_color_gradient theme_bw scale_y_discrete labs aes
#' @importFrom dplyr filter mutate
#' @importFrom stringr str_wrap
#' @importFrom forcats fct_reorder
#' @export
toppPlot <- function (toppData, category = NULL, clusters = NULL, p_val_adj="BH") {
  #check toppData object - check for columns? - check that clusters and category are in data
  test_cols = c("Category", "Name", "PValue", "GenesInTerm", "GenesInQuery", "GenesInTermInQuery")
  for (t in test_cols) {
    if (!(t %in% colnames(toppData))) {
      stop(paste0("The column ", t, " is missing from the toppData, please correct and retry."))
    }
  }
  #parse category
  if (is.null(category)){
    if (length(unique(toppData$Category)) > 1) {
      stop("Please select a ToppFun category to plot (e.g. category = 'GeneOntologyMolecularFunction').")
    }
    else if (length(unique(toppData$Category)) == 1) {
      category = unique(toppData$Category)
    }
    else {
      stop("Invalid input data - please check your toppData dataframe")
    }
  }

  #parse pvalue

  if (!(p_val_adj %in% c("BH", "Bonferroni", "BY", "none", "None"))) {
    cat("P value adjustment not found - using 'BH' by default. For no adjustment, use p_val_adj = 'none'.")
  }


  p_val_col = switch(p_val_adj,
                     "BH" = "QValueFDRBH",
                     "Bonferroni" = "QValueBonferroni",
                     "BY" = "QvalueFDRBY",
                     "none" = "PValue",
                     "None" = "PValue",
                     "QValueFDRBH")
  if (p_val_col %in% c("QValueFDRBH", "QValueBonferroni", "QvalueFDRBY")) {
    color_label = "Adj. P-value"
  } else {
    color_label = "P-value"
  }

  #parse clusters
  if (is.null(clusters)) {
    if ("Cluster" %in% colnames(toppData)) {
      clusters = unique(toppData$Cluster)
    }
  }
  if (length(clusters) > 1) {
    cat("Multiple clusters entered: function returns a list of ggplots")
    plot_list = list()
    for (c in clusters) {
      plot_list[[c]] <- toppData |>
        dplyr::filter(Cluster == c) |>
        dplyr::filter(Category == category) |>
        dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
        dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
        ggplot(mapping = aes(
          x = geneRatio,
          y = Name
        )) +
        geom_segment(aes(xend=0, yend=Name)) +
        geom_point(mapping = aes(size=GenesInTermInQuery, color=PValue)) +
        scale_color_gradient(low='red', high = 'gray') +
        theme_bw() +
        scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 20, whitespace_only = F)) +
        labs(color=color_label, size = "Genes from Query in Gene Set")
    }
    return (plot_list)
  } else if ("Cluster" %in% colnames(toppData)) {
    toppData |>
      dplyr::filter(Cluster == clusters) |>
      dplyr::filter(Category == category) |>
      dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
      dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
      ggplot(mapping = aes(
        x = geneRatio,
        y = Name
      )) +
      geom_segment(aes(xend=0, yend=Name)) +
      geom_point(mapping = aes(size=GenesInTermInQuery, color=PValue)) +
      scale_color_gradient(low='red', high = 'gray') +
      theme_bw() +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 20, whitespace_only = F)) +
      labs(color=color_label, size = "Genes from Query in Gene Set")
  } else {
    toppData |>
      dplyr::filter(Category == category) |>
      dplyr::mutate(geneRatio = GenesInTermInQuery / GenesInTerm) |>
      dplyr::mutate(Name = forcats::fct_reorder(Name, geneRatio)) |>
      ggplot(mapping = aes(
        x = geneRatio,
        y = Name
      )) +
      geom_segment(aes(xend=0, yend=Name)) +
      geom_point(mapping = aes(size=GenesInTermInQuery, color=PValue)) +
      scale_color_gradient(low='red', high = 'gray') +
      theme_bw() +
      scale_y_discrete(labels = function(x) str_wrap(x, width = 20, whitespace_only = F)) +
      labs(color=color_label, size = "Genes from Query in Gene Set")
  }

}
