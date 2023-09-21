#toppgene
library(httr)
library(rjson)

#' Get results from ToppFun
#'
#' @param markers A vector of markers or dataframe with columns as cluster labels
#' @param topp_categories A string or vector with specific toppfun categories for the query
#' @param key_type Gene name format
#' @return data.frame
#' @examples
#' CreateProjectDirectory("smith", "/Users/asmith01/projects")
#' CreateProjectDirectory("miller", "~/projects")
#' @export
ToppFun <- function(markers,
                    topp_categories,
                    key_type = "SYMBOL") {
  big_df <- data.frame()
  for (col in colnames(top_marker_df)) {
    
    if (!(col %in% c('rank', 'X'))) {
      cat('Working on cluster:', col, "\n")
      gene_list = top_marker_df[[col]]
      d <- get_topp(gene_list = gene_list,
                    topp_categories = topp_categories,
                    key_type = "SYMBOL")
      d[['cluster']] = col
      big_df <- rbind(big_df, d)
      
    }
  }
  return (big_df)
}

#' Convert genes into Entrez format
#'
#' @param genes A list of genes
#' @return a vector of genes in Entrex format
#' @examples
#' get_Entrez(genes=list_of_genes)
#' @export
get_Entrez<- function(list_of_genes){
  lookup_url = 'https://toppgene.cchmc.org/API/lookup'
  payload = toJSON(list(Symbols=list_of_genes))

  r <- POST(url = lookup_url,
            body = payload)
  new_gene_list = c()
  for (g in content(r)$Genes) {
    new_gene_list <- base::append(new_gene_list, g[['Entrez']])
  }
  return (new_gene_list)
}

get_topp <- function(gene_list,
                      key_type,
                      topp_categories,
                      max_results=10,
                      min_genes=1,
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

  data = toJSON(payload)

  #send POST request
  url = 'https://toppgene.cchmc.org/API/enrich'
  r <- POST(url = url,
            body=data)
  
  response_data <- content(r)[['Annotations']]
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

#top_markers <- read.csv('/Users/bryanwgranger/biocm/projects/ferreira/newest/main_analysis/top_markers_res07.csv')

topp_top_markers <- function(top_marker_df,
                             topp_categories) {
  big_df <- data.frame()
  for (col in colnames(top_marker_df)) {
    
    if (!(col %in% c('rank', 'X'))) {
      cat('Working on cluster:', col, "\n")
      gene_list = top_marker_df[[col]]
      d <- get_topp(gene_list = gene_list,
                    topp_categories = topp_categories,
                    key_type = "SYMBOL")
      d[['cluster']] = col
      big_df <- rbind(big_df, d)
  
    }
  }
  return (big_df)
}
#b <- topp_top_markers(top_markers, topp_categories = c("Pathway", "GeneOntologyMolecularFunction", "ToppCell"))

#' Get a vector of ToppFun categories
#'
#' @return a vector
#' @examples
#' get_ToppCats()
#' @export
get_ToppCats <- function() {
  toppCats <- c("Category","ID","Name","PValue","QValueFDRBH","QValueFDRBY","QValueBonferroni",
               "TotalGenes","GenesInTerm","GenesInQuery","GenesInTermInQuery","Source","URL")
  return (toppCats)
}