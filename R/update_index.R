## update functions for
### 1) Between-model moves for bmim
### 2) Update theta index weights

#' Update gamma_index
#'
#' Between model moves for bmim (gamma_additive, beta, gamma_nonadditive, rho, theta)
#'
#' @keywords internal
update_gamma_index <- function(params){ ##
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  x <- params$constants$x
  y_Zalpha <- params$constants$y-params$constants$z%*%params$params$alpha
  B <- params$params$B
  beta_additive <- params$params$beta_additive
  y_Bbeta_Zalpha <- y_Zalpha-B%*%beta_additive
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  theta <- params$params$theta
  xtheta <- params$params$xtheta
  prior_theta_kappa <- params$constants$prior_theta_kappa
  stepsize_theta_kappa <- params$constants$stepsize_theta_kappa
  tau2_additive <- params$params$tau2_additive
  gamma_additive <- params$params$gamma_additive
  pi_additive <- params$params$pi_additive
  tau2_nonadditive <- params$params$tau2_nonadditive
  rho_nonadditive <- params$params$rho_nonadditive
  stepsize_rho_nonadditive <- params$constants$stepsize_rho_nonadditive
  prior_rho_nonadditive <- params$constants$prior_rho_nonadditive
  gamma_nonadditive <- params$params$gamma_nonadditive
  pi_nonadditive <- params$params$pi_nonadditive
  prior_pi_nonadditive <- params$constants$prior_pi_nonadditive
  n <- params$constants$n
  p <- params$constants$p
  d <- params$constants$d
  Lm <- params$constants$Lm
  kernelfun <- params$constants$kernelfun
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank
  approxProj <- params$constants$approxProj
  SS <- params$constants$SS
  Sstar <- params$params$Sstar

  for(j in sample(1:p)){ # randomize order
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])
    beta_id <- sum(d[0:(j-1)])+(1:d[j])

    ## initialize proposals
    theta_prop <- theta
    xtheta_prop <- xtheta
    B_prop <- B
    P_prop <- P
    gamma_additive_prop <- gamma_additive
    beta_additive_prop <- beta_additive
    gamma_nonadditive_prop <- gamma_nonadditive
    rho_nonadditive_prop <- rho_nonadditive

    ## get current state
    state <- 0*(gamma_additive[j]==0 & gamma_nonadditive[j]==0)+
             1*(gamma_additive[j]==1 & gamma_nonadditive[j]==0)+
             2*(gamma_additive[j]==0 & gamma_nonadditive[j]==1)+
             3*(gamma_additive[j]==1 & gamma_nonadditive[j]==1)

    ## propose new state
    state_prop <- (state+sample(3,1))%%4 ## choose one of the other states
    gamma_additive_prop[j] <- 0+1*(state_prop==1|state_prop==3)
    gamma_nonadditive_prop[j] <- 0+1*(state_prop==2|state_prop==3)

    ## proposals and corresponding components of the MHG acceptance ratio
    logPrior <- 0
    logProp <- 0
    # proposed beta
    if(gamma_additive_prop[j]>gamma_additive[j]){
      ## Note: final term is unpenalized aka improper prior, so we draw the final term from N(0,1) separately
      beta_additive_prop[beta_id[-d[j]]] <- c(mvtnorm::rmvnorm(1,rep(0,d[j]-1),as.matrix(tau2_additive[j]*FastGP::rcppeigen_invert_matrix(Sstar[beta_id[-d[j]],beta_id[-d[j]]]))))
      beta_additive_prop[beta_id[length(beta_id)]] <- stats::rnorm(1,0,1) ## draw the final term from N(0,1) separately
      logPrior <- logPrior+  ## proposed - current
        log(pi_additive)-log(1-pi_additive) ## for gamma_additive
      ## Note: Proposal and prior cancel for beta, since we draw from the prior slab
      # logPrior <- logPrior+  ## proposed - current
      #   log(pi_additive)-log(1-pi_additive)+ ## for gamma_additive
      #   (-0.5*(d[j]-1)*log(2*pi)-0.5*(d[j]-1)*log(tau2_additive[j])+0.5*log(abs(det(as.matrix((Sstar[beta_id,beta_id])[-d[j],-d[j]]))))-0.5*as.numeric(t(beta_additive_prop[beta_id])%*%Sstar[beta_id,beta_id]%*%beta_additive_prop[beta_id]) )-0 ## for beta_additive   ## *****ADDED d[j]-1 and -length(which1) to remove unpenalized terms
      logProp <- logProp+ ## current - proposed
        0-stats::dnorm(beta_additive_prop[beta_id[length(beta_id)]],0,1,log=TRUE) ## due to final term drawn from N(0,1)
      #   0-(-0.5*(d[j]-1)*log(2*pi)-0.5*(d[j]-1)*log(tau2_additive[j])+0.5*log(abs(det(as.matrix((Sstar[beta_id,beta_id])[-d[j],-d[j]]))))-0.5*as.numeric(t(beta_additive_prop[beta_id])%*%Sstar[beta_id,beta_id]%*%beta_additive_prop[beta_id]) ) ## for beta_additive    ## *****ADDED d[j]-1 and -length(which1) to remove unpenalized terms
    }else if(gamma_additive_prop[j]<gamma_additive[j]){
      beta_additive_prop[beta_id] <- 0*beta_additive[beta_id]
      logPrior <- logPrior+  ## proposed - current
        log(1-pi_additive)-log(pi_additive) ## for gamma_additive
      ## Note: Proposal and prior cancel for beta, since we draw from the prior slab
      # logPrior <- logPrior+  ## proposed - current
      #   log(1-pi_additive)-log(pi_additive)+ ## for gamma_additive
      #   0-(-0.5*(d[j]-1)*log(2*pi)-0.5*(d[j]-1)*log(tau2_additive[j])+0.5*log(abs(det(as.matrix((Sstar[beta_id,beta_id])[-d[j],-d[j]]))))-0.5*as.numeric(t(beta_additive[beta_id])%*%Sstar[beta_id,beta_id]%*%beta_additive[beta_id]) ) ## for beta_additive    ## *****ADDED d[j]-1 and -length(which1) to remove unpenalized terms
      ## TESTING: adding in this term since acceptance ratio needs to be symmetric
      logProp <- logProp+ ## current - proposed
        stats::dnorm(beta_additive[beta_id[length(beta_id)]],0,1,log=TRUE)-0 ## due to final term drawn from N(0,1) in increasing model step above

      #   (-0.5*(d[j]-1)*log(2*pi)-0.5*(d[j]-1)*log(tau2_additive[j])+0.5*log(abs(det(as.matrix((Sstar[beta_id,beta_id])[-d[j],-d[j]]))))-0.5*as.numeric(t(beta_additive[beta_id])%*%Sstar[beta_id,beta_id]%*%beta_additive[beta_id]) )-0 ## for beta_additive     ## *****ADDED d[j]-1 and -length(which1) to remove unpenalized terms
    }
    # proposed rho
    if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){
      rho_nonadditive_prop[j] <- stats::rgamma(1,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2]) # drawing from prior
      logPrior <- logPrior+   ## proposed - current
        log(pi_nonadditive)-log(1-pi_nonadditive) ## for gamma_nonadditive
      ## Note: Proposal and prior cancel for rho, since we draw from the prior slab
      # logPrior <- logPrior+   ## proposed - current
      #   log(pi_nonadditive)-log(1-pi_nonadditive)+ ## for gamma_nonadditive
      #   dgamma(rho_nonadditive_prop[j],shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE)-0 ## for rho_nonadditive
      # logProp <- logProp+ ## current - proposed
      #   0-dgamma(rho_nonadditive_prop[j],shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE) ## for rho_nonadditive
    }else if(gamma_nonadditive_prop[j]<gamma_nonadditive[j]){
      rho_nonadditive_prop[j] <- 0
      logPrior <- logPrior+   ## proposed - current
        log(1-pi_nonadditive)-log(pi_nonadditive) ## for gamma_nonadditive
      ## Note: Proposal and prior cancel for rho, since we draw from the prior slab
      # logPrior <- logPrior+   ## proposed - current
      #   log(1-pi_nonadditive)-log(pi_nonadditive)+ ## for gamma_nonadditive
      #   0-dgamma(rho_nonadditive[j],shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE) ## for rho_nonadditive
      # logProp <- logProp+ ## current - proposed
      #   dgamma(rho_nonadditive[j],shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE)-0 ## for rho_nonadditive
    }
    # proposed theta
    if(gamma_additive[j]==0 & gamma_nonadditive[j]==0){ ## if moving from (0,0) to another state
      if(Lm[j]>1){ ## otherwise it's just 1 (same as theta)
        theta_prop[theta_id] <- c(Rfast::rvmf(1,mu=rep(1/sqrt(Lm[j]),Lm[j]),k=prior_theta_kappa))
        ## update xtheta, B, and P
        xtheta_prop[,j] <- x[[j]]%*%theta_prop[theta_id]
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
        B_prop <- newBP$B
        P_prop <- newBP$P
      }
      ## Note: Proposal and prior cancel for theta, since we draw from the prior slab
      # logPrior <- logPrior+ ## proposed - current
      #   prior_theta_kappa*sum(theta_prop[theta_id]*rep(1/sqrt(Lm[j]),Lm[j]))-0  ## for theta
      # logProp <- logProp+ ## current - proposed
      #   0-prior_theta_kappa*sum(theta_prop[theta_id]*rep(1/sqrt(Lm[j]),Lm[j]))  ## for theta
    }else if(gamma_additive_prop[j]==0 & gamma_nonadditive_prop[j]==0){ ## if moving to (0,0)
      if(Lm[j]>1){ ## otherwise it's just 1 (same as theta)
        theta_prop[theta_id] <- 0*theta[theta_id]
        ## update xtheta, B, and P
        xtheta_prop[,j] <- 0#x[[j]]%*%theta_prop[theta_id]
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
        B_prop <- newBP$B
        P_prop <- newBP$P
      }
      # ## update xtheta, B, and P
      # xtheta_prop[,j] <- 0#x[[j]]%*%theta_prop[theta_id]
      # newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
      # B_prop <- newBP$B
      # P_prop <- newBP$P
      ## Note: Proposal and prior cancel for theta, since we draw from the prior slab
      # logPrior <- logPrior+ ## proposed - current
      #   0-prior_theta_kappa*sum(theta[theta_id]*rep(1/sqrt(Lm[j]),Lm[j])) ## for theta
      # logProp <- logProp+ ## current - proposed
      #   prior_theta_kappa*sum(theta[theta_id]*rep(1/sqrt(Lm[j]),Lm[j]))-0 ## for theta
    }

    y_Bbeta_Zalpha_prop <- y_Zalpha-B_prop%*%beta_additive_prop

    ## update kernel
    K_prop <- get_K(xtheta_prop,rho_nonadditive_prop,kernelfun,invmethod,knots)
    ## update invSIG
    SIGres <- get_invSIG(P=P_prop,K=K_prop,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
    invSIG_prop <- SIGres$invSIG
    logdetSIG_prop <- SIGres$logdetSIG


    ## log of MH acceptance ratio
    logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                      as.numeric(t(y_Bbeta_Zalpha_prop)%*%(invSIG_prop)%*%y_Bbeta_Zalpha_prop-t(y_Bbeta_Zalpha)%*%(invSIG)%*%y_Bbeta_Zalpha)  )      ## proposed - current
    # logPrior <- ## proposed - current ## computed above
    # logProp <- ## current - proposed ## computed above
    logDiff <- logLik+logPrior+logProp
    logU <- log(stats::runif(1,0,1))

    if(logU < logDiff){
      theta <- theta_prop
      xtheta <- xtheta_prop
      B <- B_prop
      P <- P_prop
      K <- K_prop
      gamma_additive <- gamma_additive_prop
      beta_additive <- beta_additive_prop
      gamma_nonadditive <- gamma_nonadditive_prop
      rho_nonadditive <- rho_nonadditive_prop
      invSIG <- invSIG_prop
      logdetSIG <- logdetSIG_prop
      Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
      Dm <- SIGres$Dm
      y_Bbeta_Zalpha <- y_Bbeta_Zalpha_prop #
    }



  }




  newparams$params$theta <- theta
  newparams$params$xtheta <- xtheta
  newparams$params$B <- B
  newparams$params$P <- P
  newparams$params$K <- K
  newparams$params$gamma_additive <- gamma_additive
  newparams$params$beta_additive <- beta_additive
  newparams$params$gamma_nonadditive <- gamma_nonadditive
  newparams$params$rho_nonadditive <- rho_nonadditive
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}

update_theta <- function(params){ ## draw theta (index weights)
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  x <- params$constants$x
  y_Zalpha <- params$constants$y-params$constants$z%*%params$params$alpha
  B <- params$params$B
  beta_additive <- params$params$beta_additive
  y_Bbeta_Zalpha <- y_Zalpha-B%*%beta_additive
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  theta <- params$params$theta
  xtheta <- params$params$xtheta
  prior_theta_kappa <- params$constants$prior_theta_kappa
  stepsize_theta_kappa <- params$constants$stepsize_theta_kappa
  gamma_additive <- params$params$gamma_additive
  tau2_nonadditive <- params$params$tau2_nonadditive
  rho_nonadditive <- params$params$rho_nonadditive
  stepsize_rho_nonadditive <- params$constants$stepsize_rho_nonadditive
  prior_rho_nonadditive <- params$constants$prior_rho_nonadditive
  gamma_nonadditive <- params$params$gamma_nonadditive
  pi_nonadditive <- params$params$pi_nonadditive
  prior_pi_nonadditive <- params$constants$prior_pi_nonadditive
  n <- params$constants$n
  p <- params$constants$p
  d <- params$constants$d
  Lm <- params$constants$Lm
  kernelfun <- params$constants$kernelfun
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank
  approxProj <- params$constants$approxProj
  SS <- params$constants$SS

  ## only consider those with theta in the model, and permute ordering
  rand_scan <- (1:p)[(Lm>1) & ((gamma_additive+gamma_nonadditive)>0)]
  if(length(rand_scan)>1){
    rand_scan <- sample(rand_scan)
  }
  for(j in rand_scan){
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])

    ## initialize proposals
    theta_prop <- theta
    xtheta_prop <- xtheta

    ## proposal for theta_j
    theta_prop[theta_id] <- c(Rfast::rvmf(1,mu=theta[theta_id],k=stepsize_theta_kappa))

    ## update xtheta, B, and P
    xtheta_prop[,j] <- x[[j]]%*%theta_prop[theta_id]
    newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
    B_prop <- newBP$B
    P_prop <- newBP$P
    y_Bbeta_Zalpha_prop <- y_Zalpha-B_prop%*%beta_additive

    ## update kernel
    K_prop <- get_K(xtheta_prop,rho_nonadditive,kernelfun,invmethod,knots)
    ## update invSIG
    SIGres <- get_invSIG(P=P_prop,K=K_prop,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
    invSIG_prop <- SIGres$invSIG
    logdetSIG_prop <- SIGres$logdetSIG

    ## log of Metropolis acceptance ratio
    logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                      as.numeric(eigenQuadProd(invSIG_prop,y_Bbeta_Zalpha_prop)-eigenQuadProd(invSIG,y_Bbeta_Zalpha))  )      ## proposed - current
                      #as.numeric(t(y_Bbeta_Zalpha_prop)%*%(invSIG_prop)%*%y_Bbeta_Zalpha_prop-t(y_Bbeta_Zalpha)%*%(invSIG)%*%y_Bbeta_Zalpha)  )      ## proposed - current
    logPrior <- prior_theta_kappa*sum(theta_prop[theta_id]*rep(1/sqrt(Lm[j]),Lm[j]))-
      prior_theta_kappa*sum(theta[theta_id]*rep(1/sqrt(Lm[j]),Lm[j]))  ## proposed - current
    logDiff <- logLik+logPrior #+logProp # no proposal since using random walk metropolis (symmetry)
    logU <- log(stats::runif(1,0,1))

    if(logU < logDiff){
      theta <- theta_prop
      xtheta <- xtheta_prop
      B <- B_prop
      P <- P_prop
      K <- K_prop
      invSIG <- invSIG_prop
      logdetSIG <- logdetSIG_prop
      Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
      Dm <- SIGres$Dm
      y_Bbeta_Zalpha <- y_Bbeta_Zalpha_prop
    }


  }


  newparams$params$theta <- theta
  newparams$params$xtheta <- xtheta
  newparams$params$B <- B
  newparams$params$P <- P
  newparams$params$K <- K
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}



#' Update theta (polar)
#'
#' @keywords internal
update_theta_polar <- function(params){ ## draw theta (index weights)
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  x <- params$constants$x
  y_Zalpha <- params$constants$y-params$constants$z%*%params$params$alpha
  B <- params$params$B
  beta_additive <- params$params$beta_additive
  y_Bbeta_Zalpha <- y_Zalpha-B%*%beta_additive
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  theta <- params$params$theta
  xtheta <- params$params$xtheta
  prop_phi_a <- params$constants$prop_phi_a
  prior_theta_kappa <- params$constants$prior_theta_kappa
  stepsize_theta_kappa <- params$constants$stepsize_theta_kappa
  gamma_additive <- params$params$gamma_additive
  tau2_nonadditive <- params$params$tau2_nonadditive
  rho_nonadditive <- params$params$rho_nonadditive
  stepsize_rho_nonadditive <- params$constants$stepsize_rho_nonadditive
  prior_rho_nonadditive <- params$constants$prior_rho_nonadditive
  gamma_nonadditive <- params$params$gamma_nonadditive
  pi_nonadditive <- params$params$pi_nonadditive
  prior_pi_nonadditive <- params$constants$prior_pi_nonadditive
  n <- params$constants$n
  p <- params$constants$p
  d <- params$constants$d
  Lm <- params$constants$Lm
  kernelfun <- params$constants$kernelfun
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank
  approxProj <- params$constants$approxProj
  SS <- params$constants$SS

  ## only consider those with theta in the model, and permute ordering
  rand_scan <- (1:p)[(Lm>1) & ((gamma_additive+gamma_nonadditive)>0)]
  if(length(rand_scan)>1){
    rand_scan <- sample(rand_scan)
  }
  for(j in rand_scan){
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])

    ## polar coords
    phi_j <- get_phi(theta[theta_id])
    phibeta <- c(phi_j[1]/(pi/2),(phi_j[-1]+(pi/2))/pi) ## on beta scale for comparing proposals

    ## initialize proposals
    phi_j_prop <- phi_j
    theta_prop <- theta
    xtheta_prop <- xtheta

    ## propose new phi_j from pi*beta distribution (using previous value as mode)
    b_mode <- ((1-(phibeta))*prop_phi_a+2*(phibeta)-1)/(phibeta) ## using previous value as mode
    phibeta_prop <- stats::rbeta(length(b_mode),prop_phi_a,b_mode)
    phi_j_prop[1] <- (pi/2)*phibeta_prop[1] ## transform to appropriate scale
    phi_j_prop[-1] <- pi*phibeta_prop[-1]-(pi/2)

    theta_prop[theta_id] <- get_theta(phi_j_prop)     ## proposal for theta_j

    ## update xtheta, B, and P
    xtheta_prop[,j] <- x[[j]]%*%theta_prop[theta_id]
    newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
    B_prop <- newBP$B
    P_prop <- newBP$P
    y_Bbeta_Zalpha_prop <- y_Zalpha-B_prop%*%beta_additive

    ## update kernel
    K_prop <- get_K(xtheta_prop,rho_nonadditive,kernelfun,invmethod,knots)
    ## update invSIG
    SIGres <- get_invSIG(P=P_prop,K=K_prop,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
    invSIG_prop <- SIGres$invSIG
    logdetSIG_prop <- SIGres$logdetSIG

    ## log of Metropolis acceptance ratio
    logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                      as.numeric(eigenQuadProd(invSIG_prop,y_Bbeta_Zalpha_prop)-eigenQuadProd(invSIG,y_Bbeta_Zalpha))  )      ## proposed - current
    #as.numeric(t(y_Bbeta_Zalpha_prop)%*%(invSIG_prop)%*%y_Bbeta_Zalpha_prop-t(y_Bbeta_Zalpha)%*%(invSIG)%*%y_Bbeta_Zalpha)  )      ## proposed - current
    logPrior <- 0  ## proposed - current ## uniform prior on half sphere
    logProp <- sum(stats::dbeta(phibeta,prop_phi_a,b_mode,log=TRUE)) -
      sum(stats::dbeta(phibeta_prop,prop_phi_a,b_mode,log=TRUE)) +
      ## from change of variable
      log(abs(prod(cos(phi_j)^(Lm[j]-(1:(Lm[j]-1)))))) -
      log(abs(prod(cos(phi_j_prop)^(Lm[j]-(1:(Lm[j]-1)))))) ## current - proposed
    logDiff <- logLik+logPrior+logProp #
    logU <- log(stats::runif(1,0,1))

    if(logU < logDiff){
      theta <- theta_prop
      xtheta <- xtheta_prop
      B <- B_prop
      P <- P_prop
      K <- K_prop
      invSIG <- invSIG_prop
      logdetSIG <- logdetSIG_prop
      Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
      Dm <- SIGres$Dm
      y_Bbeta_Zalpha <- y_Bbeta_Zalpha_prop
    }


  }


  newparams$params$theta <- theta
  newparams$params$xtheta <- xtheta
  newparams$params$B <- B
  newparams$params$P <- P
  newparams$params$K <- K
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}


