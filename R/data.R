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
#' @source 10X Genomics PBMC 3k dataset. Available from 
#'   \url{https://www.10xgenomics.com/resources/datasets/}. 
#'   Analysis following Seurat PBMC tutorial: 
#'   \url{https://satijalab.org/seurat/articles/pbmc3k_tutorial.html}
#'
"pbmc.markers"

#' toppData example
#'
#' A dataframe of of sample toppData results created from the pbmc.markers dataset using the toppFun() function
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
#' @usage data("toppdata.pbmc")
#' @source \url{https://toppgene.cchmc.org}
#' @source Generated using ToppGene API (\url{https://toppgene.cchmc.org/}).
#'   Chen J, Bardes EE, Aronow BJ, Jegga AG. ToppGene Suite for gene list 
#'   enrichment analysis and candidate gene prioritization. Nucleic Acids Res. 
#'   2009;37(Web Server issue):W305-11. doi: 10.1093/nar/gkp427.
#' @source 10X Genomics PBMC 3k dataset. Available from 
#'   \url{https://www.10xgenomics.com/resources/datasets/}. 
#'   Analysis following Seurat PBMC tutorial: 
#'   \url{https://satijalab.org/seurat/articles/pbmc3k_tutorial.html}
#'
"toppdata.pbmc"
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

#' toppData example for ifnb.de
#'
#' A dataframe of of sample toppData results created from the ifnb.de dataset using the toppFun() function
#'
#' @format A dataframe with 12,227 rows and 14 columns
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
#' @usage data("toppdata.ifnb")
#' @source Generated using ToppGene API (\url{https://toppgene.cchmc.org/}).
#'   Chen J, Bardes EE, Aronow BJ, Jegga AG. ToppGene Suite for gene list 
#'   enrichment analysis and candidate gene prioritization. Nucleic Acids Res. 
#'   2009;37(Web Server issue):W305-11. doi: 10.1093/nar/gkp427.
#' @source \url{https://toppgene.cchmc.org}
#' @source Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet 
#'   single-cell RNA-sequencing using natural genetic variation. 
#'   Nat Biotechnol. 2018;36(1):89-94. doi:10.1038/nbt.4042
"toppdata.ifnb"
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
#' @source Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet 
#'   single-cell RNA-sequencing using natural genetic variation. 
#'   Nat Biotechnol. 2018;36(1):89-94. doi:10.1038/nbt.4042
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
#' @source Kang HM, Subramaniam M, Targ S, et al. Multiplexed droplet 
#'   single-cell RNA-sequencing using natural genetic variation. 
#'   Nat Biotechnol. 2018;36(1):89-94. doi:10.1038/nbt.4042
#'
"ifnb.markers.list.CD8T"

#' toppData example using the airway dataset results
#'
#' A dataframe of of sample toppData results created from the ifnb.de dataset using the toppFun() function
#'
#' @format A dataframe with 902 rows and 14 columns
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
#' @usage data("toppdata.airway")
#' @source \url{https://toppgene.cchmc.org}
#' @source Generated using ToppGene API (\url{https://toppgene.cchmc.org/}).
#'   Chen J, Bardes EE, Aronow BJ, Jegga AG. ToppGene Suite for gene list 
#'   enrichment analysis and candidate gene prioritization. Nucleic Acids Res. 
#'   2009;37(Web Server issue):W305-11. doi: 10.1093/nar/gkp427.
#' @source Himes, E. B, Jiang, X., Wagner, P., Hu, R., Wang, Q., 
#' Klanderman, B., Whitaker, M. R, Duan, Q., Lasky-Su, J., Nikolos, 
#' C., Jester, W., Johnson, M., Panettieri, A. R, Tantisira, G. K, Weiss, 
#' T. S, Lu, Q. (2014). “RNA-Seq Transcriptome Profiling Identifies CRISPLD2 
#' as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in 
#' Airway Smooth Muscle Cells.” PLoS ONE, 9(6), e99625. http://www.ncbi.nlm.nih.gov/pubmed/24926665. 
#' @source \url{https://www.bioconductor.org/packages/release/data/experiment/html/airway.html}
#'
"toppdata.airway"