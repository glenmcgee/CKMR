
#' Initialize parameters
#'
#' @keywords internal
initialize_params <- function(y, ## response
                              x, ## list of exposures/mixture components
                              z, ## confounders to be adjusted
                              prior_theta_kappa,
                              prior_pi_additive,
                              prior_tau2_additive,
                              prior_pi_nonadditive,
                              prior_rho_nonadditive,
                              prior_tau2_nonadditive,
                              prior_sigma2,
                              stepsize_theta_kappa,
                              stepsize_tau2_nonadditive,
                              stepsize_rho_nonadditive,
                              oversample,
                              prop_phi_a,
                              kernelfun,
                              invmethod,
                              rank,
                              knots, ## for GPP
                              approxProj, ## NOT USED: starting with default projection matrix anyway (default theta)
                              hierarchical,
                              polar,
                              startvals=list(), # list of starting values if given
                              basis.opts.list){
  ## constants
  n <- nrow(z) # sample size
  p <- length(x) ## (p is m now)
  pz <- ncol(z) # number of linear covariates/confounders

  ## get bases
  if(is.null(basis.opts.list)){
    basis.opts.list <- vector(mode="list",length=p)   ## if not specified, use unstructured
  }else if(!is.null(basis.opts.list$type)){
    basis.opts.list <- rep(list(basis.opts.list),p)         ## if only one option specified, use for all indices
  }
  if(length(basis.opts.list)!=p){
    basis.opts.list <- vector(mode="list",length=p)   ## if otherwise the wrong length, use unstructured
  }
  ## construct bases
  Psi <- vector(mode="list",length=p)
  for(j in (1:p)){
    basis.opts <- basis.opts.list[[j]]

    if(is.null(basis.opts)){
      Psi[[j]] <- diag(ncol(x[[j]]))                ## if not specified, use unstructured
    }else if(is.null(basis.opts$type)){
      Psi[[j]] <- diag(ncol(x[[j]]))                ## if not specified, use unstructured
    }else{
      if(toupper(basis.opts$type) %in% c("NONE","RANKED","PCA","NS","BS","FACE","GAM","MEAN","AVERAGE")){
        Psi[[j]] <- getbasis(x[[j]],basis.opts)$psi
      }else{
        stop("basis type not recognized.")
      }
    }
    x[[j]] <- x[[j]]%*%Psi[[j]] ## transform x = X Psi
  }


  ## more constants
  Lm <- sapply(x,function(a)ncol(as.matrix(a))) # index lengths



  #### params
  #### NOTE: We initialize theta twice: once for smooths; again to account for zeroes below.
  if(utils::hasName(startvals,"theta")){
    theta <- startvals$theta
  }else{
    theta <- rep(1/sqrt(Lm),times=Lm)  ## temporary in order to make smooths without zeros
  }
  xtheta <- matrix(NA,nrow=n,ncol=p)
  for(j in (1:p)){
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])
    xtheta[,j] <- as.matrix(x[[j]])%*%theta[theta_id]
  }
  ## for smooth functions
  xtheta <- data.frame(xtheta)
  # if(is.null(colnames(x))){
    colnames(xtheta) <- paste0("X",1:p)
  # }
  if(utils::hasName(startvals,"SS")){
      SS <- startvals$SS
      B <- startvals$B
      P <- startvals$P
      d <- sapply(SS,function(obj) obj$df)
  }else{
    smooth <- lapply(lapply(as.list(paste0("mgcv::s(", colnames(xtheta), ")")),str2lang),eval) ## extract list of smooths from varnames
    SS <- lapply(lapply(smooth,mgcv::smoothCon,data=xtheta,absorb.cons = TRUE),'[[',1)
    Blist <- lapply(SS,'[[','X')
    d <- sapply(Blist,ncol)
    B <- Reduce(cbind,Blist)
    Bproj <- cbind(1,B[,apply(B,2, function(x) stats::var(x)!=0)]) #
    P <- diag(1,n)-Bproj%*%MASS::ginv(t(Bproj)%*%Bproj)%*%t(Bproj) ## avoiding linear dependence issues
  }
  Sstar <- Matrix::bdiag(lapply(lapply(SS,'[[','S'),'[[',1))

  ##
  if(utils::hasName(startvals,"pi_additive")){
    pi_additive <- startvals$pi_additive
  }else{
    pi_additive <- stats::rbeta(1,prior_pi_additive[1],prior_pi_additive[2])
  }
  if(utils::hasName(startvals,"gamma_additive")){
    gamma_additive <- startvals$gamma_additive
  }else{
    gamma_additive <- stats::rbinom(p,1,pi_additive) # only one per x
  }
  if(utils::hasName(startvals,"theta")){
    theta <- startvals$theta
  }else{
    theta <- rep(gamma_additive/sqrt(Lm),times=Lm)  ## redoing theta with requisite zeros
  }
  ## updating xtheta etc. with proper gammas/thetas
  for(j in (1:p)){
    theta_id <- sum(Lm[0:(j-1)])+(1:Lm[j])
    xtheta[,j] <- as.matrix(x[[j]])%*%theta[theta_id]
    BP <- get_BP(xtheta,B,SS,j,d,n,approxProj,P)
    B <- BP$B
    P <- BP$P
  }
  ##
  if(utils::hasName(startvals,"beta_additive")){
    beta_additive <- startvals$beta_additive
  }else{
    beta_additive <- stats::rnorm(sum(d),0,1)*rep(gamma_additive,d) # should select all basis components simultaneously
  }
  if(utils::hasName(startvals,"tau2_additive")){
    tau2_additive <- startvals$tau2_additive
  }else{
    tau2_additive <- (1/stats::rgamma(p,shape=prior_tau2_additive[1],rate=prior_tau2_additive[2]))*gamma_additive
  }
  ##
  if(utils::hasName(startvals,"pi_nonadditive")){
    pi_nonadditive <- startvals$pi_nonadditive
  }else{
    pi_nonadditive <- stats::rbeta(1,prior_pi_nonadditive[1],prior_pi_nonadditive[2])
  }
  if(utils::hasName(startvals,"gamma_nonadditive")){
    gamma_nonadditive <- startvals$gamma_nonadditive
  }else{
    gamma_nonadditive <- stats::rbinom(p,1,pi_nonadditive)
  }
  if(hierarchical==TRUE){## enforce hierarchical structure
    gamma_nonadditive <- gamma_nonadditive*gamma_additive
  }

  if(utils::hasName(startvals,"rho_nonadditive")){
    rho_nonadditive <- startvals$rho_nonadditive
  }else{
    rho_nonadditive <- stats::rgamma(p,shape=prior_rho_nonadditive[1],rate=prior_rho_nonadditive[2])*gamma_nonadditive
  }
  if(utils::hasName(startvals,"tau2_nonadditive")){
    tau2_nonadditive <- startvals$tau2_nonadditive
  }else{
    tau2_nonadditive <- 1/stats::rgamma(1,shape=prior_tau2_nonadditive[1],rate=prior_tau2_nonadditive[2])
  }
  ##
  if(utils::hasName(startvals,"alpha")){
    alpha <- startvals$alpha
  }else{
    alpha <- stats::rnorm(ncol(z),0,1)
  }
  if(utils::hasName(startvals,"sigma2")){
    sigma2 <- startvals$sigma2
  }else{
    sigma2 <- 1/stats::rgamma(1,shape=prior_sigma2[1],rate=prior_sigma2[2])
  }

  ## GP covariance components
  if(invmethod=="GPP" & !is.null(rank)){
    if(is.null(knots)){
      knots  <- fields::cover.design(matrix(sapply(x,function(xx){
                                                        rep(sample(c(xx)),max(Lm))[1:(n*max(Lm))]
                                                      }),nrow=(n*max(Lm)),ncol=length(x)), nd = rank)$design
    }else{
      rank <- nrow(knots)
    }
  }else{
    knots=NULL
  }
  K <- get_K(xtheta,rho_nonadditive,kernelfun,invmethod=invmethod,knots=knots)
  SIGres <- get_invSIG(P=P,K=K,n=n,tau2_nonadditive=tau2_nonadditive,sigma2=sigma2,invmethod=invmethod,rank=rank)
  invSIG <- SIGres$invSIG
  logdetSIG <- SIGres$logdetSIG
  Um <- SIGres$Um
  Dm <- SIGres$Dm ## excludes tau

  constants <- list(n=n,
                    p=p,
                    d=d,
                    Lm=Lm,
                    pz=pz,
                    y=y,
                    x=x, ## this is the transformed (XPsi)
                    z=z,
                    knots=knots,
                    # B=B,
                    # P=P,
                    prior_theta_kappa=prior_theta_kappa,
                    prior_pi_additive=prior_pi_additive,
                    prior_tau2_additive=prior_tau2_additive,
                    prior_pi_nonadditive=prior_pi_nonadditive,
                    prior_rho_nonadditive=prior_rho_nonadditive,
                    prior_tau2_nonadditive=prior_tau2_nonadditive,
                    prior_sigma2=prior_sigma2,
                    stepsize_theta_kappa=stepsize_theta_kappa,
                    stepsize_tau2_nonadditive=stepsize_tau2_nonadditive,
                    stepsize_rho_nonadditive=stepsize_rho_nonadditive,
                    oversample=oversample,
                    prop_phi_a=prop_phi_a,
                    kernelfun=kernelfun,
                    invmethod=invmethod,
                    rank=rank,
                    approxProj=approxProj,
                    SS=SS,
                    Psi=Psi)

  params <- list(theta=theta,
                 xtheta=xtheta,
                 pi_additive=pi_additive,
                 gamma_additive=gamma_additive,
                 beta_additive=beta_additive,
                 tau2_additive=tau2_additive,
                 pi_nonadditive=pi_nonadditive,
                 gamma_nonadditive=gamma_nonadditive,
                 rho_nonadditive =rho_nonadditive ,
                 tau2_nonadditive=tau2_nonadditive,
                 alpha=alpha,
                 sigma2=sigma2,
                 K=K,
                 Um=Um,
                 Dm=Dm,
                 invSIG=invSIG,
                 logdetSIG=logdetSIG,
                 B=B,
                 P=P,
                 Sstar=Sstar)

  return(list(constants=constants,
              params=params))


}
