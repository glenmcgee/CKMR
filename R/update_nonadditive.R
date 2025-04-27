## update functions

#' Update gamma nonadditive and rho
#'
#' @keywords internal
update_gamma_nonadditive <- function(params){ ## draw both gamma and rho
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  # x <- params$constants$x
  xtheta <- params$params$xtheta
  y_Bbeta_Zalpha <- params$constants$y-
    params$params$B%*%params$params$beta-
    params$constants$z%*%params$params$alpha
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  tau2_nonadditive <- params$params$tau2_nonadditive
  rho_nonadditive <- params$params$rho_nonadditive
  stepsize_rho_nonadditive <- params$constants$stepsize_rho_nonadditive
  prior_rho_nonadditive <- params$constants$prior_rho_nonadditive
  gamma_nonadditive <- params$params$gamma_nonadditive
  pi_nonadditive <- params$params$pi_nonadditive
  prior_pi_nonadditive <- params$constants$prior_pi_nonadditive
  n <- params$constants$n
  p <- params$constants$p
  kernelfun <- params$constants$kernelfun
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank


  for(j in sample(1:p)){ # randomize order
    gamma_nonadditive_prop <- gamma_nonadditive
    rho_nonadditive_prop <- rho_nonadditive

    ## proposal for gamma[j] and rho[j]
    if(gamma_nonadditive[j]==1){ # propose 0
      gamma_nonadditive_prop[j] <- 0
      rho_nonadditive_prop[j] <- 0
    }else{
      gamma_nonadditive_prop[j] <- 1
      rho_nonadditive_prop[j] <- stats::rgamma(1,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2]) # draw from the prior (maybe we can change this)
    }

    ## update kernel
    K_prop <- get_K(xtheta,rho_nonadditive_prop,kernelfun,invmethod,knots)
    ## update invSIG
    SIGres <- get_invSIG(P=P,K=K_prop,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
    invSIG_prop <- SIGres$invSIG
    logdetSIG_prop <- SIGres$logdetSIG


    ## log of MH acceptance ratio
    logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                      as.numeric(eigenQuadProd(invSIG_prop-invSIG,y_Bbeta_Zalpha)) )
    #                    as.numeric(t(y_Bbeta_Zalpha)%*%(invSIG_prop-invSIG)%*%y_Bbeta_Zalpha) )
    ## combining prior and proposal components because of the cancellations (from drawing gamma from the prior)
    logPriorProp <- (gamma_nonadditive_prop[j])*log(pi_nonadditive/(1-pi_nonadditive)) +
                    (1-gamma_nonadditive_prop[j])*log((1-pi_nonadditive)/pi_nonadditive)
    logDiff <- logLik+logPriorProp
    logU <- log(stats::runif(1,0,1))

    if(logU < logDiff){
      gamma_nonadditive <- gamma_nonadditive_prop
      rho_nonadditive <- rho_nonadditive_prop
      K <- K_prop
      invSIG <- invSIG_prop
      logdetSIG <- logdetSIG_prop
      Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
      Dm <- SIGres$Dm
    }



  }

  newparams$params$gamma_nonadditive <- gamma_nonadditive
  newparams$params$rho_nonadditive <- rho_nonadditive
  newparams$params$K <- K
  newparams$params$invSIG <- invSIG
  newparams$params$logdetSIG <- logdetSIG
  newparams$params$Um <- Um
  newparams$params$Dm <- Dm


  return(newparams)

}

#' Refinement step for rho_nonadditive
#'
#' @keywords internal
update_rho_nonadditive <- function(params){ ## refinement step only
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  # x <- params$constants$x
  xtheta <- params$params$xtheta
  y_Bbeta_Zalpha <- params$constants$y-
    params$params$B%*%params$params$beta-
    params$constants$z%*%params$params$alpha
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  tau2_nonadditive <- params$params$tau2_nonadditive
  rho_nonadditive <- params$params$rho_nonadditive
  stepsize_tau2_nonadditive <- params$constants$stepsize_tau2_nonadditive
  stepsize_rho_nonadditive <- params$constants$stepsize_rho_nonadditive
  prior_rho_nonadditive <- params$constants$prior_rho_nonadditive
  gamma_nonadditive <- params$params$gamma_nonadditive
  n <- params$constants$n
  kernelfun <- params$constants$kernelfun
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank


  if(sum(gamma_nonadditive)>0){ ## only update if there are any non-zeroes

    if(sum(gamma_nonadditive)==1){ ## error handling since if you call (sample(5) it could return j<5)
      sample_j <- which(gamma_nonadditive!=0)
    }else{ # randomize order; only do refinement step for nonzeros
      sample_j <- sample(which(gamma_nonadditive!=0))
    }
    for(j in sample_j){
      rho_nonadditive_prop <- rho_nonadditive # initialize

      ## gamma proposal (with mean=current value, sd=stepsize)
      rho_nonadditive_prop[j] <- stats::rgamma(1,
                                        shape=rho_nonadditive[j]^2/stepsize_rho_nonadditive^2,
                                        rate=rho_nonadditive[j]/stepsize_rho_nonadditive^2)
      ## update kernel
      K_prop <- get_K(xtheta,rho_nonadditive_prop,kernelfun,invmethod,knots)
      ## update invSIG
      SIGres <- get_invSIG(P=P,K=K_prop,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
      invSIG_prop <- SIGres$invSIG
      logdetSIG_prop <- SIGres$logdetSIG

      ## log of MH acceptance ratio
      logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                        as.numeric(eigenQuadProd(invSIG_prop-invSIG,y_Bbeta_Zalpha)) )      ## proposed - current
                        #as.numeric(t(y_Bbeta_Zalpha)%*%(invSIG_prop-invSIG)%*%y_Bbeta_Zalpha) )
      # prior should just be the Gamma prior component, since this conditions on gamma_nonadditive[j]=1
      logPrior <- stats::dgamma(rho_nonadditive_prop[j],shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE)-
        stats::dgamma(rho_nonadditive[j]     ,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2],log=TRUE)                                                ## proposed - current
      logProp <- stats::dgamma(rho_nonadditive[j],     shape=rho_nonadditive_prop[j]^2/stepsize_tau2_nonadditive^2,rate=rho_nonadditive_prop[j]/stepsize_tau2_nonadditive^2)-
        stats::dgamma(rho_nonadditive_prop[j],shape=rho_nonadditive[j]^2/stepsize_tau2_nonadditive^2,     rate=rho_nonadditive[j]/stepsize_tau2_nonadditive^2)       ## current - proposed
      logDiff <- logLik+logPrior+logProp
      logU <- log(stats::runif(1,0,1))


      if(logU < logDiff){
        rho_nonadditive <- rho_nonadditive_prop
        K <- K_prop
        invSIG <- invSIG_prop
        logdetSIG <- logdetSIG_prop
        Um <- SIGres$Um ## update eigencomponents for easy use in tau2 step
        Dm <- SIGres$Dm
      }

    }

    newparams$params$rho_nonadditive <- rho_nonadditive
    newparams$params$K <- K  ## NOTE: added this back in
    newparams$params$invSIG <- invSIG
    newparams$params$logdetSIG <- logdetSIG
    newparams$params$Um <- Um
    newparams$params$Dm <- Dm
  }



  return(newparams)

}


update_pi_nonadditive <- function(params){
  newparams <- params

  p <- params$constants$p
  a_pi_rho <- params$constants$prior_pi_nonadditive[1]
  b_pi_rho <- params$constants$prior_pi_nonadditive[2]
  gamma_nonadditive <- params$params$gamma_nonadditive


  newparams$params$pi_nonadditive <- stats::rbeta(1,
                                           a_pi_rho+sum(gamma_nonadditive),
                                           b_pi_rho+p-sum(gamma_nonadditive))

  return(newparams)

}


#' Update tau2_nonadditive
#'
#' @keywords internal
update_tau2_nonadditive <- function(params){
  newparams <- params

  P <- params$params$P
  K <- params$params$K
  y_Bbeta_Zalpha <- params$constants$y-
    params$params$B%*%params$params$beta-
    params$constants$z%*%params$params$alpha
  sigma2 <- params$params$sigma2
  invSIG <- params$params$invSIG
  logdetSIG <- params$params$logdetSIG
  Um <- params$params$Um
  Dm <- params$params$Dm
  tau2_nonadditive <- params$params$tau2_nonadditive
  stepsize_tau2_nonadditive <- params$constants$stepsize_tau2_nonadditive
  prior_tau2_nonadditive <- params$constants$prior_tau2_nonadditive
  n <- params$constants$n
  invmethod <- params$constants$invmethod
  knots <- params$constants$knots
  rank <- params$constants$rank

  ## gamma proposal (with mean=current value, sd=stepsize)
  tau2_nonadditive_prop <- stats::rgamma(1,
                                  shape=tau2_nonadditive^2/stepsize_tau2_nonadditive^2,
                                  rate=tau2_nonadditive/stepsize_tau2_nonadditive^2)
  SIGres <- get_invSIG(P=P,K=K,n=n,tau2_nonadditive=tau2_nonadditive_prop,sigma2=sigma2,invmethod=invmethod,rank=rank,only_tau2=TRUE,Um=Um,Dm=Dm) ## select option to only update tau, and avoid re-approximation
  invSIG_prop <- SIGres$invSIG
  logdetSIG_prop <- SIGres$logdetSIG

  ## log of MH acceptance ratio
  logLik <- -0.5*((logdetSIG_prop-logdetSIG)+
                  as.numeric(eigenQuadProd(invSIG_prop-invSIG,y_Bbeta_Zalpha)) )                                                                                ## proposed - current
  logPrior <- stats::dgamma(1/tau2_nonadditive_prop,shape=prior_tau2_nonadditive[1],rate=prior_tau2_nonadditive[2],log=TRUE)-
    stats::dgamma(1/tau2_nonadditive     ,shape=prior_tau2_nonadditive[1],rate=prior_tau2_nonadditive[2],log=TRUE)                                        ## proposed - current
  logProp <- stats::dgamma(tau2_nonadditive     ,shape=tau2_nonadditive_prop^2/stepsize_tau2_nonadditive^2,rate=tau2_nonadditive_prop/stepsize_tau2_nonadditive^2)-
    stats::dgamma(tau2_nonadditive_prop,shape=tau2_nonadditive^2/stepsize_tau2_nonadditive^2,     rate=tau2_nonadditive/stepsize_tau2_nonadditive^2)       ## current - proposed
  logDiff <- logLik+logPrior+logProp
  logU <- log(stats::runif(1,0,1))

  if(logU < logDiff){
    newparams$params$tau2_nonadditive <- tau2_nonadditive_prop
    newparams$params$invSIG <- invSIG_prop   ## also update invSIG
    newparams$params$logdetSIG <- logdetSIG_prop
  }


  return(newparams)

}
