# Declare global variables to avoid R CMD check NOTEs
# These variables are used in ggplot2/dplyr expressions with non-standard evaluation
globalVariables(c(
    "Category", 
    "GenesInTerm", 
    "GenesInTermInQuery", 
    "Name", 
    "QValueFDRBH", 
    "geneRatio", 
    "nlog10_fdr",
    "logFC",
    "neglog_p",
    "Highlight",
    "gene"
))

# print a message warning users that data from ToppGene must be used according
# to their terms of use
.onAttach <- function(libname, pkgname) {
    msg <- "NOTE: scToppR provides data via ToppGene. Any use of this data must adhere to
  ToppGene's Terms of Use. Please visit https://toppgene.cchmc.org/navigation/termsofuse.jsp
  for more information."
    packageStartupMessage(msg)
}
