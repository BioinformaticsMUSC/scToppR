#' Get results from ToppFun
#'
#' @param markers A vector of markers or dataframe with columns as cluster labels
#' @param topp_categories A string or vector with specific toppfun categories for the query
#' @param cluster_col Column name for the groups of cells (e.g. cluster or celltype)
#' @param gene_col Column name for genes (e.g. gene or feature)
#' @param num_genes Number of genes per group to use for toppGene query (default = 20)
#' @param key_type Gene name format
#' @param p_value P-value cutoff for results
#' @param min_genes Minimum number of genes to match in a query
#' @param max_genes Maximum number of genes to match in a query
#' @param max_results Maximum number of results per cluster
#' @param correction P-value correction method ("FDR" is "BH")
#' @importFrom stringr str_glue str_c
#' @return data.frame
#' @examples
#' toppFun(markers=marker_table, topp_categories="GeneOntologyBiologicalProcess", key_type="SYMBOL")
#' @export
toppFun <- function(markers,
                    topp_categories = NULL,
                    cluster_col = "cluster",
                    gene_col = "gene",
                    num_genes = 20,
                    do_log = TRUE,
                    FDR_cutoff = 5,
                    clusters = NULL,
                    correction="FDR",
                    key_type = "SYMBOL",
                    p_value= 0.05,
                    min_genes=1,
                    max_genes=1500,
                    max_results=50
                    ) {

  #subset clusters if needed
  if (!(is.null(clusters))) {
    markers <- markers |>
      dplyr::filter(!!as.name(cluster_col) %in% clusters)
  }

  #parse input
  if ('data.frame' %in% class(markers)) {
    marker_table <- process_markers(markers=markers,
                                    cluster_col=cluster_col,
                                    gene_col=gene_col,
                                    num_genes=num_genes,
                                    FDR_cutoff=FDR_cutoff,
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

  for (col in names(marker_table)) {
    #print(col)
    if (!(col %in% c('rank', 'X'))) {
      gene_list = marker_table[[col]]
      #print(gene_list)
      if (sum(!(is.na(gene_list))) >= min_genes) {
        cat('Working on cluster:', col, '\n')
        d <- get_topp(gene_list = gene_list,
                      topp_categories = topp_categories,
                      key_type = "SYMBOL",
                      p_value=p_value,
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
  #calculate -log10(FDR)
  if (isTRUE(do_log)){
    big_df <- big_df |>
      dplyr::mutate(nlog10_fdr = -log10(QValueFDRBH))
  }
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
#' @return a vector of genes in Entrex format
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

process_markers <- function (markers, cluster_col, gene_col,
                             num_genes=20,
                             FDR_cutoff=5,
                             avg_logFC_col="avg_logFC",
                             p_val_col="p_val",
                             adj_p_val_col="p_val_adj") {


  #parse columns
  #print(colnames(markers))
  if (!(avg_logFC_col %in% colnames(markers))) {
    avg_logFC_col = stringr::str_subset(colnames(markers),
                                        pattern="avg")

  }
  if (!(adj_p_val_col %in% colnames(markers))) {
    adj_p_val_col = stringr::str_subset(colnames(markers),
                                        pattern="adj")
  }

  #avg_foldfc
  marker_table <- data.frame(row.names=seq(1, num_genes))
  for (cl in unique(markers[[cluster_col]])) {
    tdf <- markers |>
      dplyr::filter(!!as.name(cluster_col) == cl) |>
      #mutate(abs = abs(avg_log2FC)) |>
      dplyr::arrange(-!!as.name(avg_logFC_col)) |>
      dplyr::mutate(log = -log10(!!as.name(adj_p_val_col))) |>
      dplyr::filter(log > FDR_cutoff) |>
      head(num_genes) |>
      dplyr::select(!!as.name(gene_col))
    if (nrow(tdf) < num_genes) {
      new_tdf <- data.frame(row.names=1:num_genes)
      new_tdf$gene <- NA
      new_tdf[1:nrow(tdf),"gene"] = tdf[1:nrow(tdf), "gene"]
      tdf <- new_tdf
    }

    colnames(tdf) <- cl
    #rownames(tdf) <- seq(1:nrow(tdf))
    marker_table <- dplyr::bind_cols(marker_table, tdf)
  }
  #print(marker_table)
  return (marker_table)
}

get_topp <- function(gene_list,
                      key_type,
                      topp_categories,
                      max_results=10,
                      min_genes=5,
                      max_genes=1500,
                      p_value=0.05,
                      correction="FDR") {

  #assertions - to add

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
                    PValue=p_value,
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


# g = get_topp(gene_list = c("FOXP3", "IFNG"),
#              topp_categories = c("Pathway"),
#              key_type = 'SYMBOL')
#
# g = get_topp(gene_list = list(1482,4205,2626,9421,9464,6910,6722),
#              topp_categories = c("Pathway", "GeneOntologyMolecularFunction", "ToppCell"),
#              key_type = "ENTREZ")


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

