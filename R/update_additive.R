## update functions

#' Update gamma_additive
#'
#' @keywords internal
update_gamma_additive <- function(params){
  newparams <- params

  n <- params$constants$n
  d <- params$constants$d
  p <- params$constants$p
  y <- params$constants$y
  Zalpha <- params$constants$z%*%params$params$alpha
  B <- params$params$B
  Sstar <- params$params$Sstar
  gamma_additive <- params$params$gamma_additive
  beta_additive <- params$params$beta_additive
  tau2_additive <- params$params$tau2_additive
  pi_additive <- params$params$pi_additive
  invSIG <- params$params$invSIG
  invSIGmu <- invSIG%*%(y-Zalpha) ## avoiding repetitive computations

  new_gamma_additive <- gamma_additive
  for(j in sample(1:p)){ # randomize order
    gamma0 <- gamma1 <- new_gamma_additive
    gamma0[j] <- 0
    gamma1[j] <- 1
    which0 <- which(rep(gamma0,d)==1) ## only relevant columns of B
    which1 <- which(rep(gamma1,d)==1) ## only relevant columns of B
    B0 <- B[,which0]
    B1 <- B[,which1]

    V0 <- (diag(rep(gamma0/tau2_additive,d))%*%Sstar)[which0,which0]
    V1 <- (diag(rep(gamma1/tau2_additive,d))%*%Sstar)[which1,which1]
    if(length(which0)==0){ ## error handling in case there are no components left
      chol1 <- chol(as.matrix(eigenQuadProd(invSIG,B1)+V1))# as.matrix because otherwise chol gives an error about non-symmetry

      W1 <- chol2inv(chol1)

      logRratio <- 0.5*(d[j]-1)*log(tau2_additive[j])-0.5*log(det(as.matrix((Sstar[which1,which1])[-length(which1),-length(which1)])))+
         -0.5*(0-2*sum(log(diag(chol1))))+
        0.5*t(invSIGmu)%*%( - eigenQuadProd(W1,B1) )%*%invSIGmu

    }else{

      chol0 <- chol(as.matrix( eigenQuadProd(invSIG,B0)+V0))
      chol1 <- chol(as.matrix(eigenQuadProd(invSIG,B1)+V1))

      W0 <- chol2inv(chol0)
      W1 <- chol2inv(chol1)

      logRratio <- 0.5*(d[j]-1)*log(tau2_additive[j])-0.5*log(det(as.matrix((Sstar[which1[!(which1%in%which0)],which1[!(which1%in%which0)]])[-length(which1[!(which1%in%which0)]),-length(which1[!(which1%in%which0)])])))+ ## *****ADDED d[j]-1 and -length(which1) to remove unpenalized terms
        -0.5*(2*sum(log(diag(chol0)))-2*sum(log(diag(chol1))))+
        0.5*t(invSIGmu)%*%( eigenQuadProd(W0,t(B0)) - eigenQuadProd(W1,t(B1))  )%*%invSIGmu

   }

    new_gamma_additive[j] <- stats::rbinom(1,1,1/(1+((1-pi_additive)/pi_additive)*exp(logRratio)))
  }
  newparams$params$gamma_additive <- new_gamma_additive

  return(newparams)

}

#' Update beta_additive
#'
#' @keywords internal
update_beta_additive <- function(params){
  newparams <- params

  n <- params$constants$n
  d <- params$constants$d
  y <- params$constants$y
  Zalpha <- params$constants$z%*%params$params$alpha
  B <- params$params$B
  Sstar <- params$params$Sstar
  gamma_additive <- params$params$gamma_additive
  beta_additive <- params$params$beta_additive
  tau2_additive <- params$params$tau2_additive
  invSIG <- params$params$invSIG

  new_beta_additive <- 0*beta_additive

  if(sum(gamma_additive)>0){
    whichcols <- which(rep(gamma_additive,d)==1) ## only relevant columns of B
    Bgam <- B[,whichcols]
    Vgam <- (diag(rep(gamma_additive/tau2_additive,d))%*%Sstar)[whichcols,whichcols]
    V <- as.matrix(FastGP::rcppeigen_invert_matrix(as.matrix(eigenQuadProd(invSIG,Bgam)+Vgam)))
    mu <- V%*%t(Bgam)%*%invSIG%*%(y-Zalpha)

    new_beta_additive[whichcols] <- c(mvtnorm::rmvnorm(1,mu,V))

  }


  newparams$params$beta_additive <- new_beta_additive


  return(newparams)

}

#' Update pi_additive
#'
#' @keywords internal
update_pi_additive <- function(params){
  newparams <- params

  p <- params$constants$p
  a_pi <- params$constants$prior_pi_additive[1]
  b_pi <- params$constants$prior_pi_additive[2]
  gamma_additive <- params$params$gamma_additive


  newparams$params$pi_additive <- stats::rbeta(1,
                                        a_pi+sum(gamma_additive),
                                        b_pi+p-sum(gamma_additive))

  return(newparams)

}


#' Update tau2_additive
#'
#' @keywords internal
update_tau2_additive <- function(params){
  newparams <- params

  p <- params$constants$p
  d <- params$constants$d
  a_tau <- params$constants$prior_tau2_additive[1]
  b_tau <- params$constants$prior_tau2_additive[2]
  Sstar <- params$params$Sstar
  beta <- params$params$beta
  gamma_additive <- params$params$gamma_additive
  tau2_additive <- params$params$tau2_additive

  new_tau2_additive <- tau2_additive
  for(j in sample(p)){
    if(new_tau2_additive[j]!=0){ ## refinement step only
      whichx <- rep(0,p)
      whichx[j] <- 1
      whichcols <- which(rep(whichx,d)==1)
      new_tau2_additive[j] <- 1/stats::rgamma(1,
                                       shape=a_tau+gamma_additive[j]*d[j]/2,
                                       rate=b_tau+as.numeric(gamma_additive[j]*t(beta[whichcols])%*%Sstar[whichcols,whichcols]%*%beta[whichcols])/2)
    }
  }
  newparams$params$tau2_additive <- new_tau2_additive


  return(newparams)

}



