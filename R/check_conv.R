## functions for checking convergence


#' Split chains
#'
#' @keywords internal
split_chains <- function(x,nchains){
  NN <- nrow(x)/nchains
  xlist <- vector(mode = "list", length = nchains)
  for(cc in 1:nchains){
    xlist[[cc]] <- x[NN*(cc-1)+(1:NN),]
  }
  return(xlist)
}

#' Compute PSR
#'
#' @keywords internal
compute_PSR <- function(xlist){

  nchains <- length(xlist)
  NN <- nrow(xlist[[1]])   ## chain length

  ## between chain variances
  samplemeans <- Reduce('rbind',lapply(xlist,function(xx) apply(xx,2,mean)))
  B <- NN*apply(samplemeans,2,stats::var)

  ## within chain variances
  W <- apply(Reduce('rbind',lapply(xlist,function(xx) apply(xx,2,stats::var))),2,mean)

  ##
  varplus <- W*((NN-1)/NN)+B/NN ##

  return(sqrt(varplus/W))

}
