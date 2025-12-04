# print a message warning users that data from ToppGene must be used according
# to their terms of use
.onAttach <- function(libname, pkgname) {
    msg <- "NOTE: scToppR provides data via ToppGene. Any use of this data must adhere to
  ToppGene's Terms of Use. Please visit https://toppgene.cchmc.org/navigation/termsofuse.jsp
  for more information."
    packageStartupMessage(msg)
}
