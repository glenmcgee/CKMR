##

## Edited to use Rcppeigen inversion via FastGP
## Projection matrix: only use ginv (moore-penrose) if solve fails (cant use Rcppeigen here since it wont give error)


# library(Rcpp)
# library(RcppArmadillo)
# library(RcppEigen)
# library(FastGP)
# library(parallel)
#
# source("rp.R")
# source("get_K.R")
# source("get_BP.R")
# source("get_invSIG.R")
# source("initialize_params.R")
# source("update_additive.R")
# source("update_nonadditive.R")
# source("update_index.R")
# source("update_hierarchical.R")
# source("update_other.R")
# source("pred_surface.R")
# source("get_lastvals.R")
# source("get_PIPs.R")
# source("getbasis.R")
# source("get_phi.R")
# source("check_conv.R")
#
# require(mgcv)
# require(Rfast) # for rvmf
# require(MASS)


#' Fit CKMR/CMIM
#'
#' Fit collapsed kernel machine regression or multiple index model.
#'
#' @param y response
#' @param x list of exposures/mixture components
#' @param z confounders to be adjusted
#' @param niter number of iterations
#' @param nburn burn-in fraction
#' @param nthin thinning number
#' @param nchains number of chains ## not yet implemented
#' @param ncores  number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
#' @param hierarchical for BMIM, set to TRUE for hierarchically well formulated structure
#' @param polar set to TRUE for polar coordinate parameterization
#' @param prior_theta_kappa kappa hyperparameter for vMF prior (only used in bmim) (Antoniadis used 100-700)
#' @param prior_pi_additive prior for pi in spike & slab on additive component
#' @param prior_tau2_additive shape and rate for inverse gamma on spline penalty terms
#' @param prior_pi_nonadditive prior for pi in spike & slab on non-additive component
#' @param prior_rho_nonadditive shape and rate for gamma
#' @param prior_tau2_nonadditive shape and rate for inverse gamma on nonadditive GP component
#' @param prior_sigma2  shape and rate for inverse gamma on sigma^2
#' @param stepsize_theta_kappa kappa in vMF proposal for theta (bmim only)  (Antoniadis used 1000)
#' @param stepsize_tau2_nonadditive jumpsize/sd for gamma proposal on tau2_nonadditive
#' @param stepsize_rho_nonadditive jumpsize/sd for gamma proposal on rho
#' @param oversample oversample the between model moves
#' @param prop_phi_a proposal step size for phi
#' @param invmethod  "exact" for full matrix inversion, "GPP" for GP projection, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
#' @param rank  rank of approximation for lz/rp approximations of PKP (SMW) or I+PKP (DIRECT) or number of knots for GPP approach
#' @param knots  optional set of knots for GPP (overrides rank)
#' @param approxProj  set to TRUE to use the approximate projection matrix (i.e. dont update P as theta changes)
#' @param kernelfun  choice of kernel function
#' @param basis.opts.list choice of basis representations (for bkmr-DLM) #list(type="face", pve=.9)
#' @param draw_h draw h or not
#' @param centering should covariates be centered to have mean=0
#' @param scaling  should covariates be scaled to have SD=1
#' @param startvals  pass any chosen starting values
#'
#' @return Returns posterior samples
#'
#' @export
cbmim <- function(y, ## response
                  x, ## list of exposures/mixture components
                  z, ## confounders to be adjusted
                  ## MCMC specs
                  niter=10000, ## number of iterations
                  nburn=0.5*niter, ## burn-in fraction
                  nthin=10, ## thinning number
                  nchains=1, ## number of chains ## not yet implemented
                  ncores=1, ## number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
                  hierarchical=FALSE, ## for BMIM, set to TRUE for hierarchically well formulated structure
                  polar=FALSE, ## set to TRUE for polar coordinate parameterization
                  ## prior hyperparameters
                  prior_theta_kappa=1, ## kappa hyperparameter for vMF prior (only used in bmim) (Antoniadis used 100-700)
                  prior_pi_additive=c(1,1), ## prior for pi in spike & slab on additive component
                  prior_tau2_additive=c(1,1),#c(0.001,0.001), ## shape and rate for inverse gamma on spline penalty terms
                  prior_pi_nonadditive=c(1,1), ## prior for pi in spike & slab on non-additive component
                  prior_rho_nonadditive=c(1,10), ## shape and rate for gamma
                  prior_tau2_nonadditive=c(1,1),#c(0.001,0.001), ## shape and rate for inverse gamma on nonadditive GP component
                  prior_sigma2=c(1,1),#c(0.001,0.001), ## shape and rate for inverse gamma on sigma^2
                  ## MH tuning
                  stepsize_theta_kappa=100, ## kappa in vMF proposal for theta (bmim only)  (Antoniadis used 1000)
                  stepsize_tau2_nonadditive=0.1, ##jumpsize/sd for gamma proposal on tau2_nonadditive
                  stepsize_rho_nonadditive=0.1, ##jumpsize/sd for gamma proposal on rho
                  oversample=FALSE, ## oversample the between model moves
                  prop_phi_a=200,
                  ## approximating large inverse via Sherman-Morrison-Woodbury
                  invmethod=c("exact"), ## "exact" for full matrix inversion, "GPP" for GP projection, "rpSMW" for random projection +SMW, "lzSMWmax/min/both" is lanczos approximation with largest/smallest/both m eigencomponents +SMW, "lzDIRECTmax/min/both" is lanczos approximation to full matrix (I+PKP)
                  rank=100, ## rank of approximation for lz/rp approximations of PKP (SMW) or I+PKP (DIRECT) or number of knots for GPP approach
                  knots=NULL, ## optional set of knots for GPP (overrides rank)
                  approxProj=FALSE, ## set to TRUE to use the approximate projection matrix (i.e. dont update P as theta changes)
                  kernelfun="gaussian", ## choice of kernel function
                  basis.opts.list=NULL, ## choice of basis representations (for bkmr-DLM) #list(type="face", pve=.9)
                  draw_h=FALSE,
                  centering=TRUE, ## should covariates be centered to have mean=0
                  scaling=TRUE, ## should covariates be scaled to have SD=1
                  startvals=list()){ ## pass any chosen starting values

  ## add intercept since we project out a column of 1s from the kernel piece, but we are using absorb.cons=TRUE for B
  z <- cbind(1,z)

  ## define function to do analysis (called by mclapply for multiple chains)
  run_mcmc <- function(ind){
    ## initialize parameters
    params_ss <- initialize_params(y=y,x=x,z=z,
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
                                   knots=knots,
                                   approxProj=approxProj,
                                   hierarchical=hierarchical,
                                   polar=polar,
                                   startvals=startvals,
                                   basis.opts.list=basis.opts.list)

    ## define MCMC update function
    if(max(params_ss$constants$Lm)==1 & hierarchical==FALSE){ # bkmr version

      gibbs_update <- function(params){
        newparams <- params
        newparams <- update_gamma_additive(newparams) ## error due to penalization
        newparams <- update_beta_additive(newparams)
        newparams <- update_pi_additive(newparams)
        newparams <- update_tau2_additive(newparams) ## sensitive to priors?
        newparams <- update_gamma_nonadditive(newparams) ## joint step for gamma AND rho
        newparams <- update_rho_nonadditive(newparams)   ## refinement step for rho|gamma=1 to improve mixing
        newparams <- update_pi_nonadditive(newparams)
        newparams <- update_tau2_nonadditive(newparams) ## not moving at all
        newparams <- update_alpha(newparams)
        newparams <- update_sigma2(newparams) ## sensitive to priors
        return(newparams)
      }

    }else{ # bmim version

      if(hierarchical==TRUE & polar==FALSE){
        gibbs_update <- function(params){
          newparams <- params
          ### GAMMA_INDEX_HIERARCHICAL IS THE PROBLEM--> causes alpha to give warnings
          newparams <- update_gamma_index_hierarchical(newparams) ## NEW FOR HIERARCHICAL ## between model step
          newparams <- update_theta(newparams) ## refinement step for index weights|gamma=1 to improve mixing
          newparams <- update_beta_additive(newparams) ## refinement step for beta|gamma=1 to improve mixing
          newparams <- update_pi_additive(newparams)
          newparams <- update_tau2_additive(newparams) ## sensitive to priors?
          newparams <- update_rho_nonadditive(newparams)   ## refinement step for rho|gamma=1 to improve mixing
          newparams <- update_pi_nonadditive_hierarchical(newparams) ## NEW FOR HIERARCHICAL
          newparams <- update_tau2_nonadditive(newparams) ## not moving at all
          newparams <- update_alpha(newparams)
          newparams <- update_sigma2(newparams) ## sensitive to priors
          return(newparams)
        }
      }else if(hierarchical==TRUE & polar==TRUE){
        gibbs_update <- function(params){
          newparams <- params
          newparams <- update_gamma_index_hierarchical_polar(newparams) ## between model step
          newparams <- update_theta_polar(newparams) ## refinement step for index weights|gamma=1 to improve mixing
          newparams <- update_beta_additive(newparams) ## refinement step for beta|gamma=1 to improve mixing
          newparams <- update_pi_additive(newparams)
          newparams <- update_tau2_additive(newparams) ## sensitive to priors?
          newparams <- update_rho_nonadditive(newparams)   ## refinement step for rho|gamma=1 to improve mixing
          newparams <- update_pi_nonadditive_hierarchical(newparams) ##
          newparams <- update_tau2_nonadditive(newparams) ## not moving at all
          newparams <- update_alpha(newparams)
          newparams <- update_sigma2(newparams) ## sensitive to priors
          return(newparams)
        }
      }else{
        gibbs_update <- function(params){
          newparams <- params
          # newparams <- update_gamma_additive(newparams) ## error due to penalization
          # newparams <- update_gamma_nonadditive(newparams) ## joint step for gamma AND rho
          newparams <- update_gamma_index(newparams) ## between model step
          newparams <- update_theta(newparams) ## refinement step for index weights|gamma=1 to improve mixing
          newparams <- update_beta_additive(newparams) ## refinement step for beta|gamma=1 to improve mixing
          newparams <- update_pi_additive(newparams)
          newparams <- update_tau2_additive(newparams) ## sensitive to priors?
          newparams <- update_rho_nonadditive(newparams)   ## refinement step for rho|gamma=1 to improve mixing
          newparams <- update_pi_nonadditive(newparams)
          newparams <- update_tau2_nonadditive(newparams) ## not moving at all
          newparams <- update_alpha(newparams)
          newparams <- update_sigma2(newparams) ## sensitive to priors
          return(newparams)
        }
      }


    }

    ## initialize results storage
    nkeep <- round((niter-nburn)/nthin)
    keep_theta <- matrix(0,ncol=sum(params_ss$constants$Lm),nrow=nkeep)
    keep_pi_additive <- matrix(0,ncol=1,nrow=nkeep)
    keep_gamma_additive <- matrix(0,ncol=params_ss$constants$p,nrow=nkeep)
    keep_beta_additive <- matrix(0,ncol=sum(params_ss$constants$d),nrow=nkeep) ## change dimension
    keep_tau2_additive <- matrix(0,ncol=params_ss$constants$p,nrow=nkeep)
    keep_pi_nonadditive <- matrix(0,ncol=1,nrow=nkeep)
    keep_gamma_nonadditive <- matrix(0,ncol=params_ss$constants$p,nrow=nkeep)
    keep_rho_nonadditive <- matrix(0,ncol=params_ss$constants$p,nrow=nkeep)
    keep_tau2_nonadditive <- matrix(0,ncol=1,nrow=nkeep)
    keep_alpha <- matrix(0,ncol=params_ss$constants$pz,nrow=nkeep)
    keep_sigma2 <- matrix(0,ncol=1,nrow=nkeep)

    ## MCMC
    for(ss in 1:niter){

      params_ss <- gibbs_update(params_ss)

      if(ss>nburn & ss%%nthin==0){
        skeep <- (ss-nburn)/nthin
        keep_theta[skeep,] <- params_ss$params$theta
        keep_pi_additive[skeep,] <- params_ss$params$pi_additive
        keep_gamma_additive[skeep,] <- params_ss$params$gamma_additive
        keep_beta_additive[skeep,] <- params_ss$params$beta_additive
        keep_tau2_additive[skeep,] <- params_ss$params$tau2_additive
        keep_pi_nonadditive[skeep,] <- params_ss$params$pi_nonadditive
        keep_gamma_nonadditive[skeep,] <- params_ss$params$gamma_nonadditive
        keep_rho_nonadditive[skeep,] <- params_ss$params$rho_nonadditive
        keep_tau2_nonadditive [skeep,] <- params_ss$params$tau2_nonadditive
        keep_alpha[skeep,] <- params_ss$params$alpha
        keep_sigma2[skeep,] <- params_ss$params$sigma2
      }
    }

    return(list(y=y,x=x,z=z, ## note we output the RAW X (not XPsi)
                nchains=nchains,
                n=params_ss$constants$n,
                p=params_ss$constants$p,
                Lm=params_ss$constants$Lm,
                pz=params_ss$constants$pz,
                B=params_ss$params$B,
                P=params_ss$params$P,
                invmethod=invmethod,
                knots=params_ss$constants$knots,
                rank=rank,
                approxProj=approxProj,
                kernelfun=kernelfun,
                SS=params_ss$constants$SS,
                Psi=params_ss$constants$Psi,
                hierarchical=hierarchical,
                basis.opts.list=basis.opts.list,
                ##
                theta=keep_theta,
                pi_additive=keep_pi_additive,
                gamma_additive=keep_gamma_additive,
                beta_additive=keep_beta_additive,
                tau2_additive=keep_tau2_additive,
                pi_nonadditive=keep_pi_nonadditive,
                gamma_nonadditive=keep_gamma_nonadditive,
                rho_nonadditive=keep_rho_nonadditive,
                tau2_nonadditive=keep_tau2_nonadditive,
                alpha=keep_alpha,
                sigma2=keep_sigma2) )
  }

  if(nchains==1){ ## run single chain
    res <- run_mcmc()
    # res$PSR <- NULL
  }else{
    if(ncores==1){ ## run multiple chains (not parallel)
      chains <- vector(mode = "list", length = nchains)
      for(cc in 1:nchains){
        chains[[cc]] <- run_mcmc()
      }
    }else{ ## run in parallerl with nchains>1 and ncores>1
      chains <- parallel::mclapply(1:nchains,run_mcmc,mc.cores=ncores)
      # chains <- lapply(1:nchains,run_mcmc)
    }

    ## append chains
    res <- chains[[1]] ## use first for names/formatting
    collapse_names <- c("theta",            ## looping over names of variables to append
                        "pi_additive",
                        "gamma_additive",
                        "beta_additive",
                        "tau2_additive",
                        "pi_nonadditive",
                        "gamma_nonadditive",
                        "rho_nonadditive",
                        "tau2_nonadditive",
                        "alpha",
                        "sigma2")
    ## append
    res[which(names(res)%in%collapse_names)] <- lapply(collapse_names, ## looping over names of variables to append
                                                       function(v){ ## for each variable name, extract element from each chain, then collapse them via rbind
                                                         Reduce('rbind',lapply(chains,function(chn) chn[[v]]))
                                                       } )
    ## compute PSR
    res$PSR <- vector(mode = "list", length = length(collapse_names))
    res$PSR <- lapply(collapse_names,function(v){
      compute_PSR(lapply(chains,function(chn){
        chn[[v]]
      }))
    })
    names(res$PSR) <- collapse_names
  }






  return(res)

}





