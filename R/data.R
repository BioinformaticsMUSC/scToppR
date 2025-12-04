#' PBMC markers
#'
#' A dataframe of marker genes generated using the
#' FindMarkers function for each cluster from the PBMC 3k dataset
#'
#' @format A dataframe with 11,629 rows and 7 columns
#' \describe{
#'   \item{p_val}{P values}
#'   \item{avg_log2FC}{avg log 2 fc values}
#'   \item{pct.1}{percentage of cells expressing gene in group 1}
#'   \item{pct.2}{percentage of cells expressing gene in group 2}
#'   \item{p_val_adj}{adjusted p-value (FDR)}
#'   \item{cluster}{cell group name}
#'   \item{gene}{gene name}
#' }
#' @usage data("pbmc.markers")
#' @source \url{https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz}
#'
"pbmc.markers"

#' toppData example
#'
#' A dataframe of of sample toppData results
#'
#' @format A dataframe with 8,550 rows and 14 columns
#' \describe{
#'   \item{Category}{ToppGene category}
#'   \item{ID}{ToppGene Term ID}
#'   \item{Name}{ToppGene Term Name}
#'   \item{PValue}{P value}
#'   \item{QValueFDRBH}{adjusted p-value (FDR)}
#'   \item{QValueFDRBY}{adjusted p-value (BY)}
#'   \item{QValueBonferroni}{adjusted p-value (Bonferroni)}
#'   \item{TotalGenes}{Total genes in background}
#'   \item{GenesInTerm}{Genes in ToppGene Term}
#'   \item{GenesInQuery}{Genes in submitted query}
#'   \item{GenesInTermQuery}{Intersection of genes in Term and in Query}
#'   \item{Source}{ToppGene result source}
#'   \item{URL}{ToppGene associated URL}
#'   \item{Cluster}{cell group name}
#' }
#' @usage data("toppData")
#' @source \url{https://toppgene.cchmc.org}
#'
"toppData"
#' IFNB DE results
#'
#' A dataframe of differentially expressed genes generated using the
#' FindMarkers function for each cluster from the Kang 2018 IFNB dataset
#' Created using the IFNB dataset from the SeuratData package
#'
#' @format A dataframe with 92,860 rows and 7 columns
#' \describe{
#'   \item{p_val}{P values}
#'   \item{avg_log2FC}{avg log 2 fc values}
#'   \item{pct.1}{percentage of cells expressing gene in group 1}
#'   \item{pct.2}{percentage of cells expressing gene in group 2}
#'   \item{p_val_adj}{adjusted p-value (FDR)}
#'   \item{cluster}{cell group name}
#'   \item{gene}{gene name}
#' }
#' @usage data("ifnb.de")
#' @source \url{https://www.nature.com/articles/nbt.4042}
#'
"ifnb.de"

#' IFNB Marker DF
#'
#' A dataframe of 100 top markers for each class in 'seurat_annotations' column
#' using presto::wilcoxauc() and presto::top_markers()
#' Created using the IFNB dataset from the SeuratData package
#'
#' @format A dataframe with 100 rows and 14 columns
#' \describe{
#'   \item{rank}{rank of marker}
#'   \item{B}{cell group name}
#'   \item{B Activated}{cell group name}
#'   \item{CD14 Mono}{cell group name}
#'   \item{CD16 Mono}{cell group name}
#'   \item{CD4 Memory T}{cell group name}
#'   \item{CD4 Naive T}{cell group name}
#'   \item{CD8 T}{cell group name}
#'   \item{DC}{cell group name}
#'   \item{Eryth}{cell group name}
#'   \item{Mk}{cell group name}
#'   \item{CNK}{cell group name}
#'   \item{pDC}{cell group name}
#'   \item{T activated}{cell group name}
#' }
#' @usage data("ifnb.markers.df")
#' @source \url{https://www.nature.com/articles/nbt.4042}
#'
"ifnb.markers.df"

#' IFNB Marker DF
#'
#' A list of the 100 top markers for CD8 T cells in ifnb dataset
#' using presto::wilcoxauc() and presto::top_markers()
#' Created using the IFNB dataset from the SeuratData package
#'
#' @format A character vector with 100 genes
#' \describe{
#'   \item{ifnb.markers.list.CD8T}{rank of marker}
#' }
#' @usage data("ifnb.markers.list.CD8T")
#' @source \url{https://www.nature.com/articles/nbt.4042}
#'
"ifnb.markers.list.CD8T"
