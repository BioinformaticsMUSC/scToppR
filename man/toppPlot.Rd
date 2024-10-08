% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topp_plot.R
\name{toppPlot}
\alias{toppPlot}
\title{Create a dotplot from toppdata results
NEW: this accepts only 1 category}
\usage{
toppPlot(
  toppData,
  category,
  clusters = NULL,
  num_terms = 10,
  p_val_adj = "BH",
  p_val_display = "log",
  save = FALSE,
  save_dir = NULL,
  width = 5,
  height = 6,
  file_prefix = NULL,
  y_axis_text_size = 8,
  combine = FALSE,
  ncols = NULL
)
}
\arguments{
\item{toppData}{A toppData results dataframe}

\item{category}{The topp categories to plot}

\item{clusters}{The cluster(s) to plot}

\item{num_terms}{The number of terms from the toppData results to be plotted, per cluster}

\item{p_val_adj}{The P-value correction method: "BH", "Bonferroni", "BY", or "none"}

\item{p_val_display}{If "log", display the p-value in terms of -log10(p_value)}

\item{save}{Whether to save the file automatically}

\item{save_dir}{Directory to save file}

\item{width}{width of the saved file (inches)}

\item{height}{height of the saved file (inches)}

\item{file_prefix}{file prefix if saving the plot - the cluster name is also added automatically}

\item{y_axis_text_size}{Size of the Y axis text - for certain categories, it's helpful to decrease this}

\item{combine}{If TRUE and multiple clusters selected, return a patchwork object of all plots; if FALSE return list of plots}

\item{ncols}{If patchwork element returned, number of columns for subplots}
}
\value{
ggplot
}
\description{
Create a dotplot from toppdata results
NEW: this accepts only 1 category
}
\examples{
data("toppData")
toppPlot(toppData,
    category="GeneOntologyMolecularFunction",
    clusters=0,
    save=TRUE,
    file_prefix="MF_cluster0")

}
