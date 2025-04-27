## update functions for
### 1) Between-model moves for hierarchical version of bmim
### 2) Update pi_nonadditive for hierarchical version of bmim

## Between model moves for bmim (gamma_additive, beta, gamma_nonadditive, rho, theta)
#### Hierarchical version does not allow (gamma_additive,gamma_nonadditive)=(0,1)


#' Update gamma_index
#'
#' @keywords internal
update_gamma_index_hierarchical <- function(params){
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
  a_tau <- params$constants$prior_tau2_additive[1]
  b_tau <- params$constants$prior_tau2_additive[2]
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
  if(params$constants$oversample){ ## oversample larger vectors
    repeatsamp <- Lm
  }else{
    repeatsamp <- rep(1,p)
  }

  for(j in sample(1:p)){ # randomize order
    for(oversample in 1:repeatsamp[j]){ # oversample for larger vectors
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
      tau2_additive_prop <- tau2_additive

      ## get current state
      ## only allowing 3 states in hierarchical
      state <- 0*(gamma_additive[j]==0 & gamma_nonadditive[j]==0)+
        1*(gamma_additive[j]==1 & gamma_nonadditive[j]==0)+
        2*(gamma_additive[j]==1 & gamma_nonadditive[j]==1)

      ## propose new state
      ## only allowing 3 states in hierarchical
      state_prop <- (state+sample(2,1))%%3 ## choose one of the other states
      gamma_additive_prop[j] <- 0+1*(state_prop==1|state_prop==2)
      gamma_nonadditive_prop[j] <- 0+1*(state_prop==2)

      ## proposals and corresponding components of the MHG acceptance ratio
      logPrior <- 0
      logProp <- 0
      # proposed beta
      if(gamma_additive_prop[j]>gamma_additive[j]){
        tau2_additive_prop[j] <- 1/stats::rgamma(1,shape=a_tau,rate=b_tau) ## now proposing tau2_additive as well
        ## Note: final term is unpenalized aka improper prior, so we draw the final term from N(0,1) separately
        beta_additive_prop[beta_id[-d[j]]] <- c(mvtnorm::rmvnorm(1,rep(0,d[j]-1),as.matrix(tau2_additive_prop[j]*FastGP::rcppeigen_invert_matrix(as.matrix(Sstar[beta_id[-d[j]],beta_id[-d[j]]])))))
        beta_additive_prop[beta_id[length(beta_id)]] <- stats::rnorm(1,0,1) ## draw the final term from N(0,1) separately
        logPrior <- logPrior+  ## proposed - current
          log(pi_additive)-log(1-pi_additive) ## for gamma_additive
        ## Note: Proposal and prior cancel for beta, since we draw from the prior slab
          logProp <- logProp+ ## current - proposed
          0-stats::dnorm(beta_additive_prop[beta_id[length(beta_id)]],0,1,log=TRUE) ## due to final term drawn from N(0,1)
        ## gamma_nonadditive (hierarchical)
        if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){##0-->1
          logPrior <- logPrior+  ## proposed - current
            log(pi_nonadditive) ## for gamma_nonadditive
        }else{##0-->0
          logPrior <- logPrior+  ## proposed - current
            log(1-pi_nonadditive)
        }
      }else if(gamma_additive_prop[j]<gamma_additive[j]){
        beta_additive_prop[beta_id] <- 0*beta_additive[beta_id]
        tau2_additive_prop[j] <- 0
        logPrior <- logPrior+  ## proposed - current
          log(1-pi_additive)-log(pi_additive) ## for gamma_additive
        ## Note: Proposal and prior cancel for beta, since we draw from the prior slab
        logProp <- logProp+ ## current - proposed
          stats::dnorm(beta_additive[beta_id[length(beta_id)]],0,1,log=TRUE)-0 ## due to final term drawn from N(0,1) in increasing model step above
        ## gamma_nonadditive (hierarchical)
        if(gamma_nonadditive_prop[j]<gamma_nonadditive[j]){##1-->0
          logPrior <- logPrior+  ## proposed - current
            -log(pi_nonadditive) ## for gamma_nonadditive
        }else{##0-->0
          logPrior <- logPrior+  ## proposed - current
            -log(1-pi_nonadditive) ## for gamma_nonadditive
        }
      }
      ## gamma_nonadditive if gamma_additive remained at 1
      if(gamma_additive_prop[j]==1 & gamma_additive[j]==1){
        if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){##0-->1
          logPrior <- logPrior+  ## proposed - current
            log(pi_nonadditive)-log(1-pi_nonadditive)
        }else{##1-->0
          logPrior <- logPrior+  ## proposed - current
            -log(pi_nonadditive)+log(1-pi_nonadditive)
        }
      }
      # proposed rho
      if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){
        rho_nonadditive_prop[j] <- stats::rgamma(1,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2]) # drawing from prior
        ## NOTE: already added prior for gamma_nonadditive above in hierarchical version
      }else if(gamma_nonadditive_prop[j]<gamma_nonadditive[j]){
        rho_nonadditive_prop[j] <- 0
        ## NOTE: already added prior for gamma_nonadditive above in hierarchical version
      }
      # proposed theta
      if(gamma_additive[j]==0 & gamma_nonadditive[j]==0){ ## if moving from (0,0) to another state
        ## update xtheta, B, and P
        if(Lm[j]>1){ ## otherwise it's just 1 (same as theta)
          theta_prop[theta_id] <- c(Rfast::rvmf(1,mu=rep(1/sqrt(Lm[j]),Lm[j]),k=prior_theta_kappa))
          xtheta_prop[,j] <- x[[j]]%*%theta_prop[theta_id]
        }else{
          theta_prop[theta_id] <- 1
          xtheta_prop[,j] <- x[[j]]
        }
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)## can speed this up later by avoiding new P computation
        B_prop <- newBP$B
        P_prop <- newBP$P

      }else if(gamma_additive_prop[j]==0 & gamma_nonadditive_prop[j]==0){ ## if moving to (0,0)
        theta_prop[theta_id] <- 0*theta[theta_id]
        ## update xtheta, B, and P
        xtheta_prop[,j] <- 0
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
        B_prop <- newBP$B
        P_prop <- newBP$P

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
                        as.numeric(eigenQuadProd(invSIG_prop,y_Bbeta_Zalpha_prop)-eigenQuadProd(invSIG,y_Bbeta_Zalpha))  )      ## proposed - current
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
        tau2_additive <- tau2_additive_prop
        invSIG <- invSIG_prop
        logdetSIG <- logdetSIG_prop
        Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
        Dm <- SIGres$Dm
        y_Bbeta_Zalpha <- y_Bbeta_Zalpha_prop #
      }


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
  newparams$params$tau2_additive <- tau2_additive
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}


#' Update gamma_index_hierarchical (Polar)
#'
#' @keywords internal
update_gamma_index_hierarchical_polar <- function(params){
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
  a_tau <- params$constants$prior_tau2_additive[1]
  b_tau <- params$constants$prior_tau2_additive[2]
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
  if(params$constants$oversample){ ## oversample larger vectors
    repeatsamp <- Lm
  }else{
    repeatsamp <- rep(1,p)
  }

  for(j in sample(1:p)){ # randomize order
    for(oversample in 1:repeatsamp[j]){ # oversample for larger vectors
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
      tau2_additive_prop <- tau2_additive

      ## get current state
      ## only allowing 3 states in hierarchical
      state <- 0*(gamma_additive[j]==0 & gamma_nonadditive[j]==0)+
        1*(gamma_additive[j]==1 & gamma_nonadditive[j]==0)+
        2*(gamma_additive[j]==1 & gamma_nonadditive[j]==1)

      ## propose new state
      ## only allowing 3 states in hierarchical
      state_prop <- (state+sample(2,1))%%3 ## choose one of the other states
      gamma_additive_prop[j] <- 0+1*(state_prop==1|state_prop==2)
      gamma_nonadditive_prop[j] <- 0+1*(state_prop==2)

      ## proposals and corresponding components of the MHG acceptance ratio
      logPrior <- 0
      logProp <- 0
      # proposed beta
      if(gamma_additive_prop[j]>gamma_additive[j]){
        tau2_additive_prop[j] <- 1/stats::rgamma(1,shape=a_tau,rate=b_tau) ## EDIT: now proposing tau2_additive as well
        ## Note: final term is unpenalized aka improper prior, so we draw the final term from N(0,1) separately
        beta_additive_prop[beta_id[-d[j]]] <- c(mvtnorm::rmvnorm(1,rep(0,d[j]-1),as.matrix(tau2_additive_prop[j]*FastGP::rcppeigen_invert_matrix(as.matrix(Sstar[beta_id[-d[j]],beta_id[-d[j]]])))))
        beta_additive_prop[beta_id[length(beta_id)]] <- stats::rnorm(1,0,1) ## draw the final term from N(0,1) separately
        logPrior <- logPrior+  ## proposed - current
          log(pi_additive)-log(1-pi_additive) ## for gamma_additive
        logProp <- logProp+ ## current - proposed
          0-stats::dnorm(beta_additive_prop[beta_id[length(beta_id)]],0,1,log=TRUE) ## due to final term drawn from N(0,1)
        ## gamma_nonadditive (hierarchical)
        if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){##0-->1
          logPrior <- logPrior+  ## proposed - current
            log(pi_nonadditive) ## for gamma_nonadditive
        }else{##0-->0
          logPrior <- logPrior+  ## proposed - current
            log(1-pi_nonadditive)
        }
      }else if(gamma_additive_prop[j]<gamma_additive[j]){
        beta_additive_prop[beta_id] <- 0*beta_additive[beta_id]
        tau2_additive_prop[j] <- 0## EDIT: now proposing tau2_additive as well
        logPrior <- logPrior+  ## proposed - current
          log(1-pi_additive)-log(pi_additive) ## for gamma_additive
        logProp <- logProp+ ## current - proposed
          stats::dnorm(beta_additive[beta_id[length(beta_id)]],0,1,log=TRUE)-0 ## due to final term drawn from N(0,1) in increasing model step above

        ## gamma_nonadditive (hierarchical)
        if(gamma_nonadditive_prop[j]<gamma_nonadditive[j]){##1-->0
          logPrior <- logPrior+  ## proposed - current
            -log(pi_nonadditive) ## for gamma_nonadditive
        }else{##0-->0
          logPrior <- logPrior+  ## proposed - current
            -log(1-pi_nonadditive)
        }
      }
      ## gamma_nonadditive if gamma_additive remained at 1
      if(gamma_additive_prop[j]==1 & gamma_additive[j]==1){
        if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){##0-->1
          logPrior <- logPrior+  ## proposed - current
            log(pi_nonadditive)-log(1-pi_nonadditive)
        }else{##1-->0
          logPrior <- logPrior+  ## proposed - current
            -log(pi_nonadditive)+log(1-pi_nonadditive)
        }
      }
      # proposed rho
      if(gamma_nonadditive_prop[j]>gamma_nonadditive[j]){
        rho_nonadditive_prop[j] <- stats::rgamma(1,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2]) # drawing from prior
      }else if(gamma_nonadditive_prop[j]<gamma_nonadditive[j]){
        rho_nonadditive_prop[j] <- 0
      }
      # proposed theta
      if(gamma_additive[j]==0 & gamma_nonadditive[j]==0){ ## if moving from (0,0) to another state
        ## update xtheta, B, and P
        if(Lm[j]>1){ ## otherwise it's just 1 (same as theta)
          theta_prop[theta_id] <- c(Rfast::rvmf(1,mu=rep(1,Lm[j]),k=10e-200)) ## draw from uniform prior
          theta_prop[theta_id[1]] <- abs(theta_prop[theta_id[1]]) ## first element positive
          xtheta_prop[,j] <- x[[j]]%*%theta_prop[theta_id]
        }else{
          theta_prop[theta_id] <- 1
          xtheta_prop[,j] <- x[[j]]
        }
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)## can speed this up later by avoiding new P computation
        B_prop <- newBP$B
        P_prop <- newBP$P

      }else if(gamma_additive_prop[j]==0 & gamma_nonadditive_prop[j]==0){ ## if moving to (0,0)
        theta_prop[theta_id] <- 0*theta[theta_id]
        ## update xtheta, B, and P
        xtheta_prop[,j] <- 0
        newBP <- get_BP(xtheta_prop,B,SS,j,d,n,approxProj,P)
        B_prop <- newBP$B
        P_prop <- newBP$P

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
                        as.numeric(eigenQuadProd(invSIG_prop,y_Bbeta_Zalpha_prop)-eigenQuadProd(invSIG,y_Bbeta_Zalpha))  )      ## proposed - current
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
        tau2_additive <- tau2_additive_prop
        invSIG <- invSIG_prop
        logdetSIG <- logdetSIG_prop
        Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
        Dm <- SIGres$Dm
        y_Bbeta_Zalpha <- y_Bbeta_Zalpha_prop #
      }


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
  newparams$params$tau2_additive <- tau2_additive
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}


#' Update pi_nonadditive
#'
#' @keywords internal
update_pi_nonadditive_hierarchical <- function(params){
  newparams <- params

  p <- params$constants$p
  a_pi_rho <- params$constants$prior_pi_nonadditive[1]
  b_pi_rho <- params$constants$prior_pi_nonadditive[2]
  gamma_additive <- params$params$gamma_additive
  gamma_nonadditive <- params$params$gamma_nonadditive


  newparams$params$pi_nonadditive <- stats::rbeta(1,
                                           a_pi_rho+sum(gamma_additive*gamma_nonadditive),
                                           b_pi_rho+sum(gamma_additive*(1-gamma_nonadditive)))

  return(newparams)

}



