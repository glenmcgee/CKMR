

#' Get Kernel matrix
#'
#' @keywords internal
get_K <- function(xtheta,rho_nonadditive,kernelfun,invmethod,knots){

  if(kernelfun=="gaussian"){
    if(invmethod=="GPP"){
      Xtheta1rho <- sweep(knots, 2, sqrt(rho_nonadditive), "*")
      Xtheta0rho <- sweep(xtheta, 2, sqrt(rho_nonadditive), "*")
      K11 <- exp(-fields::rdist(Xtheta1rho, Xtheta1rho)^2 )
      K10 <- exp(-fields::rdist(Xtheta1rho, Xtheta0rho)^2 )
      return(list(K=NULL,K11=K11,K10=K10))
    }else{
      Xthetarho <- sweep(xtheta, 2, sqrt(rho_nonadditive), `*`)
      K <- exp(-fields::rdist(Xthetarho, Xthetarho)^2)
      return(list(K=K,K11=NULL,K10=NULL))
    }
  }

}
