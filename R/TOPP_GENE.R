#' Get results from ToppFun
#'
#' @param markers A vector of markers or dataframe with columns as cluster labels
#' @param topp_categories A string or vector with specific toppfun categories for the query
#' @param cluster_col Column name for the groups of cells (e.g. cluster or celltype)
#' @param gene_col Column name for genes (e.g. gene or feature)
#' @param num_genes Number of genes per group to use for toppGene query
#' @param key_type Gene name format
#' @param pval_cutoff (adjusted) P-value cutoff for filtering differentially expressed genes
#' @param fc_filter Include "ALL" genes, or only "UPREG" or "DOWNREG" for each cluster
#' @param fc_cutoff Avg log fold change cutoff for filtering differentially expressed genes
#' @param clusters Which clusters to include in toppGene query
#' @param min_genes Minimum number of genes to match in a query
#' @param max_genes Maximum number of genes to match in a query
#' @param max_results Maximum number of results per cluster
#' @param correction P-value correction method ("FDR" is "BH")
#' @importFrom stringr str_glue str_c str_subset
#' @importFrom dplyr filter arrange select
#' @return data.frame
#' @examples
#' toppFun(markers=marker_table, topp_categories="GeneOntologyBiologicalProcess", key_type="SYMBOL")
#' @export
toppFun <- function(markers,
                    topp_categories = NULL,
                    cluster_col = "cluster",
                    gene_col = "gene",
                    num_genes = NULL,
                    pval_cutoff = 0.5,
                    fc_cutoff = 0,
                    fc_filter = "ALL",
                    genes_submit_cutoff = 1000,
                    clusters = NULL,
                    correction="FDR",
                    key_type = "SYMBOL",
                    min_genes=1,
                    max_genes=1500,
                    max_results=50
                    ) {

  #subset clusters if needed
  if (!(is.null(clusters))) {
    markers <- markers |>
      dplyr::filter(!!as.name(cluster_col) %in% clusters)
  }

  #parse fc_filter
  if (!(fc_filter %in% c("ALL", "UPREG", "DOWNREG"))){
    stop("please select one of c('ALL', 'UPREG', 'DOWNREG') for fc_filter")
  }

  #parse input
  if ('data.frame' %in% class(markers)) {
    marker_list <- process_markers(markers=markers,
                                    cluster_col=cluster_col,
                                    gene_col=gene_col,
                                    num_genes=num_genes,
                                    pval_cutoff=pval_cutoff,
                                   fc_cutoff = fc_cutoff,
                                   genes_submit_cutoff=genes_submit_cutoff,
                                    fc_filter=fc_filter,
                                    avg_logFC_col="avg_logFC",
                                    p_val_col="p_val",
                                    adj_p_val_col="p_val_adj")
  } else {
    stop("data format not recognized")
  }

  #parse correction method
  if (!(correction %in% c("none", "FDR", "Bonferroni"))) {
    stop("invalid P-value correction method, please select either 'none', 'FDR', or 'Bonferroni'")
  }

  #parse categories
  if (is.null(topp_categories)) {
    topp_categories = get_ToppCats()
  }

  big_df <- data.frame()
  missing_clusters = c()

  for (col in names(marker_list)) {
    #print(col)
    if (!(col %in% c('rank', 'X'))) {
      gene_list = marker_list[[col]]
      #print(gene_list)
      #print(gene_list)
      if (sum(!(is.na(gene_list))) >= min_genes) {
        cat('Working on cluster:', col, '\n')
        d <- get_topp(gene_list = gene_list,
                      topp_categories = topp_categories,
                      key_type = "SYMBOL",
                      pval_cutoff=pval_cutoff,
                      min_genes=min_genes,
                      max_genes=max_genes,
                      max_results=max_results,
                      correction=correction)
        #print(d)
        if (nrow(d) == 0){
          missing_clusters = append(missing_clusters, col)
        } else if (length(names(markers)) == 1) {
          big_df <- rbind(big_df, d)
        } else {
          d[['Cluster']] = col
          big_df <- rbind(big_df, d)
        }
      } else {
        missing_clusters = append(missing_clusters, col)
      }
    }
  }
  #print(big_df)

  #report any missing clusters
  if (length(missing_clusters) > 0) {
    write(stringr::str_glue("WARNING: no results found for clusters {stringr::str_c(missing_clusters, collapse=', ')}"),
          stderr())
  }
  return (big_df)
}

#' Convert genes into Entrez format
#'
#' @param genes A list of genes
#' @return a vector of genes in Entrez format
#' @importFrom rjson toJSON
#' @importFrom httr POST content
#' @examples
#' get_Entrez(genes=c("IFNG", "FOXP3"))
#' @export
get_Entrez<- function(genes){
  lookup_url = 'https://toppgene.cchmc.org/API/lookup'
  payload = rjson::toJSON(list(Symbols=genes))

  r <- httr::POST(url = lookup_url,
            body = payload)
  new_gene_list = c()

  for (g in httr::content(r)$Genes) {
    new_gene_list <- base::append(new_gene_list, g[['Entrez']])
  }
  return (new_gene_list)
}

##### PROCESS MARKER INPUTS
process_markers <- function (markers, cluster_col, gene_col,
                             num_genes=1000,
                             pval_cutoff=0.5,
                             fc_cutoff=0,
                             fc_filter="ALL",
                             genes_submit_cutoff=1000,
                             avg_logFC_col="avg_logFC",
                             p_val_col="p_val",
                             adj_p_val_col="p_val_adj") {


  #parse columns - tries to find columns if differ from the default
  if (!(avg_logFC_col %in% colnames(markers))) {
    avg_logFC_col = stringr::str_subset(colnames(markers),
                                        pattern="avg")

  }
  if (!(adj_p_val_col %in% colnames(markers))) {
    adj_p_val_col = stringr::str_subset(colnames(markers),
                                        pattern="adj")
  }

  #Make a list of lists - each sub list has all of the genes (up to num_genes if specified) of the filtered data
  marker_list = list()

  for (cl in unique(markers[[cluster_col]])) {
    all_cl_markers <- markers |>
      dplyr::filter(!!as.name(cluster_col) == cl) |>
      dplyr::filter(!!as.name(adj_p_val_col) < pval_cutoff)
    if (fc_filter == "ALL"){
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(abs(!!as.name(avg_logFC_col)) > fc_cutoff) |>
        dplyr::arrange(-abs(!!as.name(avg_logFC_col))) |>
        dplyr::select(!!as.name(gene_col))
    } else if (fc_filter == "UPREG") {
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(!!as.name(avg_logFC_col) > fc_cutoff) |>
        dplyr::arrange(-!!as.name(avg_logFC_col)) |>
        dplyr::select(!!as.name(gene_col))
    } else if (fc_filter == "DOWNREG") {
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(!!as.name(avg_logFC_col) < fc_cutoff) |>
        dplyr::arrange(!!as.name(avg_logFC_col)) |>
        dplyr::select(!!as.name(gene_col))
    }
    if (!(is.null(num_genes))) {
      if (length(all_cl_markers[[gene_col]]) > num_genes) {
        marker_list[[cl]] <- all_cl_markers[1:num_genes, gene_col] |> unlist() |> as.character()
      } else {
        marker_list[[cl]] <- all_cl_markers[[gene_col]] |> unlist() |> as.character()
      }
    } else {
      marker_list[[cl]] <- all_cl_markers[[gene_col]] |> unlist() |> as.character()
    }
    }
  return (marker_list)
}

get_topp <- function(gene_list,
                      key_type,
                      topp_categories,
                      max_results=10,
                      min_genes=5,
                      max_genes=1500,
                      pval_cutoff=0.05,
                      correction="FDR") {

  #assertions - to add
  #print(gene_list)
  #convert gene names if necessary
  if (key_type != 'ENTREZ') {
    new_gene_list = get_Entrez(gene_list)
  } else{
    new_gene_list = gene_list
  }

  #create payload for POST request
  payload = list()
  payload[['Genes']] = new_gene_list
  category_list = list()
  for (i in 1:length(topp_categories)) {
    cat_dict = list(Type=topp_categories[i],
                    PValue=pval_cutoff,
                    MinGenes=min_genes,
                    MaxGenes=max_genes,
                    MaxResults=max_results,
                    Correction=correction)
    category_list[[i]] <- cat_dict
  }

  payload[['Categories']] = category_list

  data = rjson::toJSON(payload)

  #send POST request
  url = 'https://toppgene.cchmc.org/API/enrich'
  r <- httr::POST(url = url,
            body=data)

  response_data <- httr::content(r)[['Annotations']]
  return_df <- NULL
  keepers <- c("Category","ID","Name","PValue","QValueFDRBH","QValueFDRBY","QValueBonferroni",
               "TotalGenes","GenesInTerm","GenesInQuery","GenesInTermInQuery","Source","URL")
  for (i in 1:length(response_data)) {
    if (is.null(return_df)) {
      return_df <- data.frame(response_data[[i]][keepers])
    } else {
      return_df <- rbind(return_df,response_data[[i]][keepers])
    }
  }
  return (return_df)
}


#' Get a vector of ToppFun categories
#'
#' @return a vector
#' @examples
#' get_ToppCats()
#' @export
get_ToppCats <- function() {
  toppCats <- c("GeneOntologyMolecularFunction",
                "GeneOntologyBiologicalProcess",
                "GeneOntologyCellularComponent",
                "HumanPheno",
                "MousePheno",
                "Domain",
                "Pathway",
                "Pubmed",
                "Interaction",
                "Cytoband",
                "TFBS",
                "GeneFamily",
                "Coexpression",
                "CoexpressionAtlas",
                "ToppCell",
                "Computational",
                "MicroRNA",
                "Drug",
                "Disease")
  return (toppCats)
}

#' Save toppData results (optionally) split by celltype/cluster
#'
#' @param toppData Results from toppFun as a dataframe
#' @param filename filename prefix for each split file
#' @param split Boolean, whether to split the dataframe by celltype/cluster
#' @param format Saved file format, one of c("xlsx", "csv", "tsv")
#' @importFrom openxlsx write.xlsx
#' @importFrom stringr str_glue
#' @importFrom dplyr filter
#' @examples
#' toppSave(toppData, filename="toppFun_results", split = TRUE, format = "xlsx")
#' @export
toppSave <- function (toppData,
                      filename = NULL,
                      split = TRUE,
                      format = "xlsx") {
  if (isTRUE(split)) {
    for (gr in unique(toppData$Cluster)) {
      tmp_toppData <- toppData |>
        dplyr::filter(Cluster == gr)

      if (format == 'xlsx') {

        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData_{gr}.xlsx")
        } else {
          this_file = stringr::str_glue("{filename}_{gr}.xlsx")
        }

        openxlsx::write.xlsx(tmp_toppData,
                             file = this_file,
                             colNames = TRUE,
                             rowNames = FALSE,
                             borders = "columns",
                             sheetName="toppData",
                             overwrite=T)
        cat("Saving file:", this_file, "\n")

      } else if (format == "csv") {
        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData_{gr}.csv")
        } else {
          this_file = stringr::str_glue("{filename}_{gr}.csv")
        }

        write.csv(tmp_toppData,
                  file = this_file,
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      } else if (format == "tsv") {
      if (is.null(filename)) {
        this_file = stringr::str_glue("toppData_{gr}.tsv")
      } else {
        this_file = stringr::str_glue("{filename}_{gr}.tsv")
      }

        write.table(tmp_toppData,
                  file = this_file,
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      }
    }
  } else {
      if (format == 'xlsx') {

        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData.xlsx")
        } else {
          this_file = stringr::str_glue("{filename}.xlsx")
        }

        openxlsx::write.xlsx(tmp_toppData,
                             file = this_file,
                             colNames = TRUE,
                             rowNames = FALSE,
                             borders = "columns",
                             sheetName="toppData",
                             overwrite=T)
        cat("Saving file:", this_file, "\n")

      } else if (format == "csv") {
        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData.csv")
        } else {
          this_file = stringr::str_glue("{filename}.csv")
        }

        write.csv(tmp_toppData,
                  file = this_file,
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      } else if (format == "tsv") {
      if (is.null(filename)) {
        this_file = stringr::str_glue("toppData.tsv")
      } else {
        this_file = stringr::str_glue("{filename}.tsv")
      }

      write.table(tmp_toppData,
                  file = this_file,
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
      cat("Saving file:", this_file, "\n")
    }

  }

}
