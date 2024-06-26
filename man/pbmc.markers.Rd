% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pbmc.markers}
\alias{pbmc.markers}
\title{PBMC markers}
\format{
A dataframe with 11,629 rows and 7 columns
\describe{
\item{p_val}{P values}
\item{avg_log2FC}{avg log 2 fc values}
\item{pct.1}{percentage of cells expressing gene in group 1}
\item{pct.2}{percentage of cells expressing gene in group 2}
\item{p_val_adj}{adjusted p-value (FDR)}
\item{cluster}{cell group name}
\item{gene}{gene name}
}
}
\source{
\url{https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
}
\usage{
data("pbmc.markers")
}
\description{
A dataframe of marker genes generated using the
FindMarkers function for each cluster from the PBMC 3k dataset
}
\keyword{datasets}
