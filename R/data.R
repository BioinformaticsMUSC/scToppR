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
