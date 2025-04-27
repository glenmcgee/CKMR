
#' Compute theta from phi
#'
#' @keywords internal
get_theta <- function(phi){
  return(c(sin(phi),1) * cumprod(c(1,cos(phi))))
}

#' Compute phi from theta
#'
#' @keywords internal
get_phi <- function(theta){

  phi <- c(asin(theta[1]))
  if(length(theta)>2){
    for(jj in 2:(length(theta)-1)){
      phi <- c(phi,
                min(max(asin(theta[jj]/(prod(cos(phi)))),-1),1) ## in case of rounding error
               )
    }
  }
  if(sum(is.na(phi))){
    print(theta)
    print(phi)
    print(cos(phi))
  }
  return(phi)
}



