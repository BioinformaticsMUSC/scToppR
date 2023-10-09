#' Get results from ToppFun
#'
#' @param markers A vector of markers or dataframe with columns as cluster labels
#' @param topp_categories A string or vector with specific toppfun categories for the query
#' @param key_type Gene name format
#' @param p_value P-value cutoff for results
#' @param min_genes Minimum number of genes to match in a query
#' @param max_genes Maximum number of genes to match in a query
#' @param max_results Maximum number of results per cluster
#' @param correction P-value correction method ("FDR" is "BH")
#' @importFrom stringr str_glue str_c
#' @return data.frame
#' @examples
#' toppFun(markers=c("IFNG", "FOXP3"), topp_categories="GeneOntologyBiologicalProcess", key_type="SYMBOL")
#' @export
toppFun <- function(markers,
                    topp_categories,
                    key_type = "SYMBOL",
                    p_value= 0.05,
                    min_genes=1,
                    max_genes=1500,
                    max_results=10,
                    correction="FDR") {
  big_df <- data.frame()
  missing_clusters = c()
  for (col in colnames(markers)) {


    if (!(col %in% c('rank', 'X'))) {
      cat('Working on cluster:', col, "\n")
      gene_list = markers[[col]]
      d <- get_topp(gene_list = gene_list,
                    topp_categories = topp_categories,
                    key_type = "SYMBOL",
                    p_value=p_value,
                    min_genes=min_genes,
                    max_genes=max_genes,
                    max_results=max_results,
                    correction=correction)
      if (nrow(d) == 0){
        missing_clusters = append(missing_clusters, col)
      } else {
        d[['Cluster']] = col
        big_df <- rbind(big_df, d)
      }
    }

  }
  #report any missing clusters
  if (length(missing_clusters) > 0) {
    write(stringr::str_glue("WARNING: no results found for clusters {stringr::str_c(missing_clusters, collapse=',')}"),
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
