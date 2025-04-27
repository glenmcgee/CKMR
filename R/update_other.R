## update functions

#' Update alpha
#'
#' @keywords internal
update_alpha <- function(params){
  newparams <- params

  n <- params$constants$n
  pz <- params$constants$pz
  y <- params$constants$y
  Bbeta <- params$params$B%*%params$params$beta
  Z <- params$constants$z
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG


  V <- FastGP::rcppeigen_invert_matrix(diag(1,pz)+t(Z)%*%invSIG%*%Z)#solve(diag(1,pz)+t(Z)%*%invSIG%*%Z)
  mu <- V%*%t(Z)%*%invSIG%*%(y-Bbeta)


  newparams$params$alpha <- c(mvtnorm::rmvnorm(1,mu,V))


  return(newparams)

}

#' Update sigma^2
#'
#' @keywords internal
update_sigma2 <- function(params){
  newparams <- params

  a_sig <- params$constants$prior_sigma2[1]
  b_sig <- params$constants$prior_sigma2[2]
  n <- params$constants$n
  y <- params$constants$y
  Bbeta <- params$params$B%*%params$params$beta
  Zalpha <- params$constants$z%*%params$params$alpha
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG

  newsigma2 <- 1/stats::rgamma(1,
                        shape=a_sig+0.5*n,
                        # rate=b_sig+0.5*t(y-Bbeta-Zalpha)%*%solve(sigma2*invSIG,(y-Bbeta-Zalpha))) ## note multiplying by sigma2 because invSIG is divided by sigma2
                        rate=b_sig+0.5*t(y-Bbeta-Zalpha)%*%(sigma2*invSIG)%*%(y-Bbeta-Zalpha)) ## note multiplying by sigma2 because invSIG is divided by sigma2

  newparams$params$sigma2 <- newsigma2
  newparams$params$invSIG <- sigma2*invSIG/newsigma2 ## update invSIG as well because it includes the sigma scaling factor
  newparams$params$logdetSIG <- logdetSIG-n*log(sigma2)+n*log(newsigma2)

  return(newparams)
}

