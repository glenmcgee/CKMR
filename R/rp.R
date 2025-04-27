## source Rcpp function from https://github.com/jwpark88/projMCML
## m doesnt do anything
## need to truncate at m
## set k=2*m
## set alpha=1
## reports d, but need to square them posthoc
# library(Rcpp)
# library(RcppArmadillo)
# library(RcppEigen)

#' RP approximation
#'
#' @keywords internal
Rcpp::sourceCpp("./src/rp.cpp")
