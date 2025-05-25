
#' Predict surface
#'
#' Predict regression surface, componentwise by default, holding other exposures at 0
#'
#' @param fit Fitted model
#' @param Xnew New exposure values for manual/arbitrary prediction
#' @param gridlen Length of exposure grid
#' @param includeInt  Include intercept in prediction
#'
#' @return Returns list of predictions
#'
#' @export
pred_surface <- function(fit,
                         Xnew=NULL,  ## for manual/arbitrary prediction
                         gridlen=21,
                         includeInt=TRUE){ ## include intercept in prediction

  ##
  n <- fit$n
  p <- fit$p
  y <- fit$y
  x <- fit$x ## original x (not transformed)
  z <- fit$z
  S <- length(fit$sigma2)
  kernelfun <- fit$kernelfun
  invmethod <- fit$invmethod
  knots <- fit$knots
  rank <- fit$rank
  SS <- fit$SS
  theta <- fit$theta
  beta_additive <- fit$beta_additive
  alpha <- fit$alpha
  rho_nonadditive <- fit$rho_nonadditive
  tau2_nonadditive <- fit$tau2_nonadditive
  sigma2 <- fit$sigma2
  Psi=fit$Psi

  # grid of new exposures
  Lm <- sapply(x,function(a)ncol(as.matrix(a))) # original Lms (temporary)
  if(is.null(Xnew)){
    componentwise=TRUE
    Xnew <- vector(mode = "list", length = p)
    for(jj in 1:p){
      Xnew[[jj]] <- matrix(0,byrow=TRUE,ncol=Lm[jj],nrow=sum(Lm)*gridlen) ## add new rows of 0 at the mean
      for(ll in 1:Lm[jj]){
        Xnew[[jj]][gridlen*(sum(Lm[(1:jj)-1])+ll-1)+1:gridlen,ll] <- seq(-2,2,length=gridlen)
      }
    }
  }else{
    componentwise=FALSE
  }

  ## apply bases
  for(jj in 1:p){
    x[[jj]] <- x[[jj]]%*%Psi[[jj]]
    Xnew[[jj]] <- Xnew[[jj]]%*%Psi[[jj]]
  }
  Lm <- fit$Lm ## now get proper Lms for transformed matrices (x=XPsi)


  # initialize xtheta
  nnew <- nrow(as.matrix(Xnew[[1]]))
  nfull <- n+nnew
  Xthetafull <- matrix(NA,nrow=nfull,ncol=p)
  Xthetafull <- data.frame(Xthetafull)
  colnames(Xthetafull) <- paste0("X",1:p)


  pred_mat <- matrix(NA,nrow=S,ncol=nnew)
  for(ss in 1:S){ ## loop over draws

    for(jj in (1:p)){
      theta_id <- sum(Lm[0:(jj-1)])+(1:Lm[jj])
      Xthetafull[,jj] <- rbind(as.matrix(x[[jj]]),
                          as.matrix(Xnew[[jj]]))%*%theta[ss,theta_id]
    }

    Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod,knots)
    if(is.null(Kfull$K)){
      Kfull$K <- t(Kfull$K10)%*%FastGP::rcppeigen_invert_matrix(Kfull$K11)%*%Kfull$K10
      if(anyNA(Kfull$K)){
        Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod="exact",knots)
      }
    }


    Bfull <- Reduce(cbind,lapply(SS,mgcv::PredictMat,data=data.frame(Xthetafull)))
    Bfullproj <- cbind(1,Bfull[,apply(Bfull,2, function(x) stats::var(x)!=0)]) ## remove columns that are constant due to gamma_additive=0 (hence non invertible)
    Bfull <- cbind(1,Bfull) ## add intercept for projection only
    B <- Bfull[1:n,-1]
    Bnew <- Bfull[(n+1):nfull,-1]

    try(invBB <- solve(eigenMapMatMult(t(Bfullproj),Bfullproj)), ## catching linear dependence errors
        invBB <- MASS::ginv(eigenMapMatMult(t(Bfullproj),Bfullproj)))
    Pfull <- diag(1,nfull)-eigenMapMatMult(Bfullproj,invBB%*%t(Bfullproj))

    PKPtfull <- eigenQuadProd(Kfull$K,t(Pfull))#Pfull%*%Kfull%*%t(Pfull)
    PKPt <- PKPtfull[1:n,1:n]
    PKnnPt <- PKPtfull[(n+1):nfull,(n+1):nfull]
    PKnoPt <- PKPtfull[(n+1):nfull,1:n]

    #
    W <- FastGP::rcppeigen_invert_matrix(diag(1,n)+tau2_nonadditive[ss]*PKPt)#solve(diag(1,n)+tau2_nonadditive[ss]*PKPt)

    hmean <- tau2_nonadditive[ss]*eigenMapMatMult(eigenMapMatMult(PKnoPt,(W)),(y-B%*%beta_additive[ss,]-z%*%alpha[ss,])) ## multiplying by sigma2 since sigma2 gets included in invSIG
    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*eigenQuadProd((W),t(PKnoPt)))#    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*PKnoPt%*%(W)%*%t(PKnoPt))
    try(
      draw_nonadditive <- mgcv::rmvn(1,c(hmean),(hcov+t(hcov))/2)
    )

    draw_additive <- Bnew%*%beta_additive[ss,]

    pred_mat[ss,] <- c(draw_additive)+draw_nonadditive

    if(includeInt==TRUE){
      pred_mat[ss,] <- pred_mat[ss,]+alpha[ss,1]
    }
  }

  if(componentwise==TRUE){
    ## make into list for easier extraction
    pred_list <- vector(mode = "list", length = p)
    for(jj in 1:sum(Lm)){ ## loop over exposures
      pred_list[[jj]] <-  pred_mat[,gridlen*(jj-1)+(1:gridlen)]#
    }
    return(pred_list)
  }else{
    return(pred_mat)
  }



}


#' Predict surface
#'
#' Predict regression surface, indexwise, holding other indices at 0
#'
#' @param fit Fitted model
#' @param gridlen Length of exposure grid
#' @param includeInt  Include intercept in prediction
#'
#' @return Returns list of predictions
#'
#' @export
pred_surface_indexwise <- function(fit,
                                   # Xnew=NULL,
                                   gridlen=21,
                                   includeInt=TRUE){ ## include intercept in prediction

  ##
  n <- fit$n
  p <- fit$p
  Lm <- fit$Lm
  y <- fit$y
  x <- fit$x
  z <- fit$z
  S <- length(fit$sigma2)
  kernelfun <- fit$kernelfun
  invmethod <- fit$invmethod
  knots <- fit$knots
  rank <- fit$rank
  SS <- fit$SS
  theta <- fit$theta
  beta_additive <- fit$beta_additive
  alpha <- fit$alpha
  rho_nonadditive <- fit$rho_nonadditive
  tau2_nonadditive <- fit$tau2_nonadditive
  sigma2 <- fit$sigma2
  Psi <- fit$Psi

  # posterior mean of theta
  meantheta <- apply(theta,2,mean)

  ## apply bases (to original x and it also gets applied to Xnew (which is built from x))
  for(jj in 1:p){
    x[[jj]] <- x[[jj]]%*%Psi[[jj]]
  }



  # new xtheta
  Xthetanew <- matrix(0,nrow=p*gridlen,ncol=p)
  Xthetanew <- data.frame(Xthetanew)
  colnames(Xthetanew) <- paste0("X",1:p)
  nfull <- n+nrow(Xthetanew)
  for(jj in (1:p)){
    theta_id <- sum(Lm[0:(jj-1)])+(1:Lm[jj])
    thetabar <- meantheta[theta_id]/sqrt(sum(meantheta[theta_id]^2))
    thetabar[is.na(thetabar)] <- 1/sqrt(length(thetabar[is.na(thetabar)]^2)) ## if everything is zero
    Xthetanew[gridlen*(jj-1)+1:gridlen,jj] <- seq(stats::quantile(x[[jj]]%*%(thetabar),0.05),
                            stats::quantile(x[[jj]]%*%(thetabar),0.95),
                            length=gridlen)
  }
  # initialize xtheta
  Xthetafull <- matrix(NA,nrow=nfull,ncol=p)
  Xthetafull <- data.frame(Xthetafull)
  colnames(Xthetafull) <- paste0("X",1:p)



  pred_mat <- matrix(NA,nrow=S,ncol=gridlen*p)
  for(ss in 1:S){ ## loop over draws

    for(jj in (1:p)){
      theta_id <- sum(Lm[0:(jj-1)])+(1:Lm[jj])
      Xthetafull[,jj] <- c(as.matrix(x[[jj]])%*%theta[ss,theta_id], ## observed xtheta, using theta
                          Xthetanew[,jj]) ## add the new grid, which only uses meantheta
    }


    Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod,knots)
    if(is.null(Kfull$K)){
      Kfull$K <- t(Kfull$K10)%*%FastGP::rcppeigen_invert_matrix(Kfull$K11)%*%Kfull$K10
      if(anyNA(Kfull$K)){
        Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod="exact",knots)
      }
    }

    Bfull <- Reduce(cbind,lapply(SS,mgcv::PredictMat,data=data.frame(Xthetafull)))
    Bfull <- cbind(1,Bfull) ## add intercept for projection only
    B <- Bfull[1:n,-1]
    Bnew <- Bfull[(n+1):nfull,-1]
    Bfullproj <- cbind(1, ## add intercept for projection only
                       Bfull[,apply(Bfull,2, function(x) stats::var(x)!=0)]) ## remove columns that are constant due to gamma_additive=0 (hence non invertible)

    try(invBB <- solve(eigenMapMatMult(t(Bfullproj),Bfullproj)), ## catching linear dependence errors
        invBB <- MASS::ginv(eigenMapMatMult(t(Bfullproj),Bfullproj)))
    Pfull <- diag(1,nfull)-eigenMapMatMult(Bfullproj,invBB%*%t(Bfullproj))
    PKPtfull <- eigenQuadProd(Kfull$K,t(Pfull))
    PKPt <- PKPtfull[1:n,1:n]
    PKnnPt <- PKPtfull[(n+1):nfull,(n+1):nfull]
    PKnoPt <- PKPtfull[(n+1):nfull,1:n]

    #
    W <- FastGP::rcppeigen_invert_matrix(diag(1,n)+tau2_nonadditive[ss]*PKPt)#solve(diag(1,n)+tau2_nonadditive[ss]*PKPt)

    hmean <- tau2_nonadditive[ss]*eigenMapMatMult(eigenMapMatMult(PKnoPt,(W)),(y-B%*%beta_additive[ss,]-z%*%alpha[ss,]))#hmean <- tau2_nonadditive[ss]*PKnoPt%*%(W)%*%(y-B%*%beta_additive[ss,]-z%*%alpha[ss,]) ## multiplying by sigma2 since sigma2 gets included in invSIG
    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*eigenQuadProd((W),t(PKnoPt)))#    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*PKnoPt%*%(W)%*%t(PKnoPt))
    try(
      draw_nonadditive <- mgcv::rmvn(1,c(hmean),hcov)
    )
    if(!exists("draw_nonadditive")){
      draw_nonadditive <- mgcv::rmvn(1,c(hmean),
                               (hcov+t(hcov))/2) ## if it is non-symmetric
    }

    draw_additive <- Bnew%*%beta_additive[ss,]

    pred_mat[ss,] <- draw_additive+draw_nonadditive

    if(includeInt==TRUE){
      pred_mat[ss,] <- pred_mat[ss,]+alpha[ss,1]
    }
  }

  ## make into list for easier extraction
  pred_list <- vector(mode = "list", length = p)
  for(jj in 1:p){ ## loop over exposures
    pred_list[[jj]] <-  pred_mat[,gridlen*(jj-1)+(1:gridlen)]#
  }

  ##
  return(pred_list)


}



#' Summarize prediction intervals for single matrix ofresults to be called by summarize_pred().
#' @keywords internal
summarize_pred_mat <- function(pred,
                               pct=0.95){## for CIs

  res <- cbind(apply(pred,2,mean),
               apply(pred,2,stats::sd),
               apply(pred,2,function(x) stats::quantile(x,(1-pct)/2)),
               apply(pred,2,function(x) stats::quantile(x,pct+(1-pct)/2)))
  colnames(res) <- c("mean","se","lci","uci")

  ## centered
  ones <- rep(1,ncol(pred))
  H <- ones%*%solve(t(ones)%*%ones)%*%t(ones) ## projection onto column of 1s
  I_H <- diag(1,nrow=nrow(H))-H
  pred_cent <- t(I_H%*%t(pred))
  centered <- cbind(apply(pred_cent,2,mean),
               apply(pred_cent,2,stats::sd),
               apply(pred_cent,2,function(x) stats::quantile(x,(1-pct)/2)),
               apply(pred_cent,2,function(x) stats::quantile(x,pct+(1-pct)/2)))
  colnames(centered) <- c("mean","se","lci","uci")

  ## associations/contrast
  if(ncol(pred)%%2==0){ ## if no median, use adjacent quantile
    pred_assoc <- pred-matrix(pred[,ncol(pred)/2],ncol=ncol(pred),nrow=nrow(pred))
  }else{ ## compared to median median
    pred_assoc <- pred-matrix(pred[,(ncol(pred)+1)/2],ncol=ncol(pred),nrow=nrow(pred))
  }
  assoc <- cbind(apply(pred_assoc,2,mean),
               apply(pred_assoc,2,stats::sd),
               apply(pred_assoc,2,function(x) stats::quantile(x,(1-pct)/2)),
               apply(pred_assoc,2,function(x) stats::quantile(x,pct+(1-pct)/2)))
  colnames(assoc) <- c("mean","se","lci","uci")

  return(list(fitted=data.frame(res),
              centered=data.frame(centered),
              assoc=data.frame(assoc)))
}



#' Sumarize prediction intervals
#'
#' Summarize prediction intervals results of a call to pred_surface
#'
#' @param pred Predictions
#' @param pct Pct level for CIs
#' @param assoc Estimate contrast (T/F)
#' @param centered Center estimates (T/F) only used if assoc=FALSE
#'
#' @return Returns list of predictions
#'
#' @export
summarize_pred <- function(pred,
                           pct=0.95,## for CIs
                           assoc=TRUE,
                           centered=TRUE){ ## only used if assoc=FALSE
  if(assoc==TRUE){
    if(is.list(pred)){
      if(is.list(pred[[1]])){
        return(lapply(pred,function(x) lapply(x,function(y) lapply(y,function(z) summarize_pred_mat(z,pct)$assoc))))
      }else{
        return(lapply(pred,function(x) summarize_pred_mat(x,pct)$assoc))
      }

    }else{
      return(summarize_pred_mat(pred,pct)$assoc)
    }
  }else if(centered==TRUE){
    if(is.list(pred)){
      if(is.list(pred[[1]])){
        return(lapply(pred,function(x) lapply(x,function(y) lapply(y,function(z) summarize_pred_mat(z,pct)$centered))))
      }else{
        return(lapply(pred,function(x) summarize_pred_mat(x,pct)$centered))
      }
    }else{
      return(summarize_pred_mat(pred,pct)$centered)
    }
  }else{
    if(is.list(pred)){
      if(is.list(pred[[1]])){
        return(lapply(pred,function(x) lapply(x,function(y) lapply(y,function(z) summarize_pred_mat(z,pct)$fitted))))
      }else{
        return(lapply(pred,function(x) summarize_pred_mat(x,pct)$fitted))
      }
    }else{
      return(summarize_pred_mat(pred,pct)$fitted)
    }
  }


}


#' Summarize PSR for a call to pred_surface
#' @keywords internal
summarize_pred_PSR <- function(pred,nchains){

  if(is.list(pred)){
    if(is.list(pred[[1]])){
      if(is.list(pred[[1]][[1]])){
        return(lapply(pred,function(x) lapply(x,function(y) lapply(y,function(z) compute_PSR(split_chains(z,nchains))))))
      }else{
        return(lapply(pred,function(x) lapply(x,function(y) compute_PSR(split_chains(y,nchains)))))
      }
    }else{
      return(lapply(pred,function(x) compute_PSR(split_chains(x,nchains))))
    }
  }else{
    return(compute_PSR(split_chains(pred,nchains)))
  }

}





#' Indexwise predictions with interactions
#'
#' @param fit Fitted model
#' @param gridlen Length of exposure grid
#' @param includeInt  Include intercept in prediction
#' @param whichids Which ids to plot (overrides restrict option)
#' @param restrict restrict to non-null associations (only for interactions)
#'
#' @return Returns list of predictions
#'
#' @export
pred_surface_indexwise_interac <- function(fit,
                                   # Xnew=NULL,
                                   gridlen=21,
                                   includeInt=TRUE,## include intercept in prediction
                                   whichids=NULL,
                                   restrict=TRUE ){## restrict to non-nulls (only for interactions)


  ##
  n <- fit$n
  p <- fit$p
  Lm <- fit$Lm
  y <- fit$y
  x <- fit$x
  z <- fit$z
  S <- length(fit$sigma2)
  kernelfun <- fit$kernelfun
  invmethod <- fit$invmethod
  knots <- fit$knots
  rank <- fit$rank
  SS <- fit$SS
  theta <- fit$theta
  beta_additive <- fit$beta_additive
  alpha <- fit$alpha
  rho_nonadditive <- fit$rho_nonadditive
  tau2_nonadditive <- fit$tau2_nonadditive
  sigma2 <- fit$sigma2
  Psi <- fit$Psi
  if(!is.null(whichids)){
    whichindex <- whichids
  }else if(restrict==TRUE){
    whichindex <- which(apply(fit$gamma_nonadditive,2,mean)>0.5)
  }else{
    whichindex <- 1:p
  }


  # posterior mean of theta
  meantheta <- apply(theta,2,mean)

  ## apply bases (to original x and it also gets applied to Xnew (which is built from x))
  for(jj in 1:p){
    x[[jj]] <- x[[jj]]%*%Psi[[jj]]
  }



  # new xtheta
  numcompar <- length(whichindex)*(length(whichindex)-1)*3
  Xthetanew <- matrix(0,nrow=numcompar*gridlen,ncol=p)
  Xthetanew <- data.frame(Xthetanew)
  colnames(Xthetanew) <- paste0("X",1:p)
  nfull <- n+nrow(Xthetanew)
  for(jj_iter in (1:length(whichindex))){
      jj <- whichindex[jj_iter]
      theta_id <- sum(Lm[0:(jj-1)])+(1:Lm[jj])
      thetabar <- meantheta[theta_id]/sqrt(sum(meantheta[theta_id]^2))
      thetabar[is.na(thetabar)] <- 1/sqrt(length(thetabar[is.na(thetabar)]^2)) ## if everything is zero
      Xthetanew[((length(whichindex)-1)*3*gridlen)*(jj_iter-1)+1:((length(whichindex)-1)*3*gridlen),jj] <- rep(seq(stats::quantile(x[[jj]]%*%(thetabar),0.05),
                                                    stats::quantile(x[[jj]]%*%(thetabar),0.95),
                                                    length=gridlen),(length(whichindex)-1)*3)
      for(kk_iter in (1:(length(whichindex)-1))){
          kk <- (whichindex[-jj_iter])[kk_iter]
          theta_id <- sum(Lm[0:(kk-1)])+(1:Lm[kk])
          thetabar <- meantheta[theta_id]/sqrt(sum(meantheta[theta_id]^2))
          thetabar[is.na(thetabar)] <- 1/sqrt(length(thetabar[is.na(thetabar)]^2)) ## if everything is zero
          Xthetanew[(((length(whichindex)-1)*3*gridlen)*(jj_iter-1)+1:((length(whichindex)-1)*3*gridlen))[(3*gridlen)*(kk_iter-1)+1:(3*gridlen)]  , ## of the relevant grid ox Xjs, which correspond to Xks
                    kk] <- rep(c(stats::quantile(x[[kk]]%*%(thetabar),0.1),
                                   stats::quantile(x[[kk]]%*%(thetabar),0.5),
                                    stats::quantile(x[[kk]]%*%(thetabar),0.9)),each=gridlen)




      }


  }
  # initialize xtheta
  Xthetafull <- matrix(NA,nrow=nfull,ncol=p)
  Xthetafull <- data.frame(Xthetafull)
  colnames(Xthetafull) <- paste0("X",1:p)



  pred_mat <- matrix(NA,nrow=S,ncol=nrow(Xthetanew))
  for(ss in 1:S){ ## loop over draws

    for(jj in (1:p)){
      theta_id <- sum(Lm[0:(jj-1)])+(1:Lm[jj])
      Xthetafull[,jj] <- c(as.matrix(x[[jj]])%*%theta[ss,theta_id], ## observed xtheta, using theta
                           Xthetanew[,jj]) ## add the new grid, which only uses meantheta
    }

    Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod,knots)
    if(is.null(Kfull$K)){
      Kfull$K <- t(Kfull$K10)%*%FastGP::rcppeigen_invert_matrix(Kfull$K11)%*%Kfull$K10
      if(anyNA(Kfull$K)){
        Kfull <- get_K(Xthetafull,rho_nonadditive[ss,],kernelfun,invmethod="exact",knots)
      }
    }

    Bfull <- Reduce(cbind,lapply(SS,mgcv::PredictMat,data=data.frame(Xthetafull)))
    Bfull <- cbind(1,Bfull) ## add intercept for projection only
    B <- Bfull[1:n,-1]
    Bnew <- Bfull[(n+1):nfull,-1]
    Bfullproj <- cbind(1, ## add intercept for projection onl
                       Bfull[,apply(Bfull,2, function(x) stats::var(x)!=0)]) ## remove columns that are constant due to gamma_additive=0 (hence non invertible)

    try(invBB <- solve(eigenMapMatMult(t(Bfullproj),Bfullproj)), ## catching linear dependence errors
        invBB <- MASS::ginv(eigenMapMatMult(t(Bfullproj),Bfullproj)))
    Pfull <- diag(1,nfull)-eigenMapMatMult(Bfullproj,invBB%*%t(Bfullproj))
    PKPtfull <- eigenQuadProd(Kfull$K,t(Pfull))
    PKPt <- PKPtfull[1:n,1:n]
    PKnnPt <- PKPtfull[(n+1):nfull,(n+1):nfull]
    PKnoPt <- PKPtfull[(n+1):nfull,1:n]

    # # not using the approximation because
    W <- FastGP::rcppeigen_invert_matrix(diag(1,n)+tau2_nonadditive[ss]*PKPt)#solve(diag(1,n)+tau2_nonadditive[ss]*PKPt)


    hmean <- tau2_nonadditive[ss]*eigenMapMatMult(eigenMapMatMult(PKnoPt,(W)),(y-B%*%beta_additive[ss,]-z%*%alpha[ss,]))#hmean <- tau2_nonadditive[ss]*PKnoPt%*%(W)%*%(y-B%*%beta_additive[ss,]-z%*%alpha[ss,]) ## multiplying by sigma2 since sigma2 gets included in invSIG
    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*eigenQuadProd((W),t(PKnoPt)))#    hcov <- sigma2[ss]*tau2_nonadditive[ss]*(PKnnPt-tau2_nonadditive[ss]*PKnoPt%*%(W)%*%t(PKnoPt))
    try(
      draw_nonadditive <- mgcv::rmvn(1,c(hmean),hcov)
    )
    if(!exists("draw_nonadditive")){
      draw_nonadditive <- mgcv::rmvn(1,c(hmean),
                               (hcov+t(hcov))/2) ## if it is non-symmetric
    }

    draw_additive <- Bnew%*%beta_additive[ss,]

    pred_mat[ss,] <- draw_additive+draw_nonadditive

    if(includeInt==TRUE){
      pred_mat[ss,] <- pred_mat[ss,]+alpha[ss,1]
    }
  }

  ## make into list for easier extraction
  pred_list <- vector(mode = "list", length = length(whichindex))
  overall_iter <- 0
  for(jj_iter in (1:length(whichindex))){
    jj <- whichindex[jj_iter]
      pred_list[[jj_iter]] <-  vector(mode = "list", length = (length(whichindex)-1))
      for(kk_iter in (1:(length(whichindex)-1))){
        kk <- (whichindex[-jj_iter])[kk_iter]
        pred_list[[jj_iter]][[kk_iter]] <-  vector(mode = "list", length = 3)
            for(aa in 1:3){
              overall_iter <- overall_iter+1
              pred_list[[jj_iter]][[kk_iter]][[aa]] <- pred_mat[,gridlen*(overall_iter-1)+(1:gridlen)]#
            }



        }


  }

  ##
  return(pred_list)


}


