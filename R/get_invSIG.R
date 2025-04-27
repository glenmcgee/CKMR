
## function to update invSIG and logdetSIG
## only using invmethod="exact" right now

#' Get Sigma inverse and log of det(Sigma)
#'
#' @keywords internal
get_invSIG <- function(P,K,n=nrow(P),tau2_nonadditive,sigma2,invmethod="exact",rank=NULL,
                       only_tau2=FALSE,Um=NULL,Dm=NULL){ ## need to pass Um and Dm if using the only_tau2 option to avoid re-approximating PKP

  if(invmethod=="exact"){

    SIG <- sigma2*(diag(1,n)+tau2_nonadditive*eigenQuadProd(K$K,P))#
    invSIG <- FastGP::rcppeigen_invert_matrix(SIG)
    cholSIG <- FastGP::rcppeigen_get_chol(SIG)
    logdetSIG <- 2*as.numeric(sum(log((diag(cholSIG)))))
    Um <- Dm <- NULL

  }else if(invmethod=="GPP"){ ## GP projection + Sherman Morrison Woodbury as in Savitsky et al and Bobb et al 2018

    n0 <- ncol(K$K10)
    n1 <- nrow(K$K10)
    Q <- K$K11 + diag(0.001, n1, n1)
    R <- Q + tau2_nonadditive*tcrossprod(eigenMapMatMult(K$K10,P))
    cholQ <- FastGP::rcppeigen_get_chol(Q)
    cholR <- FastGP::rcppeigen_get_chol(R)
    Rinv <- chol2inv(t(cholR)) ## need the transpose because rcpp eigen gives lower triangular version
    invSIG <- (diag(1, n0, n0) - tau2_nonadditive*eigenQuadProd(Rinv,eigenMapMatMult(K$K10,P)))/sigma2 ##
    logdetSIG <- n0*log(sigma2)-(2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR))))
    Um <- Dm <- NULL

  }else if(invmethod=="rpSMW"){ ## random projection + Sherman Morrison Woodbury as in Guan et al and Banerjee et al

    if(only_tau2==TRUE & !is.null(Um) & !is.null(Dm)){ ## if only updating tau2, we can skip the approximation and just manipulate tau2
      Um <- Um
      Dm <- Dm
    }else{
      rpm <- rp(n,rank*2,1,eigenQuadProd(K$K,P),rank)#rp(n,rank*2,1,t(P)%*%K$K%*%P,rank)
      Dm <- (c(rpm$d)[1:rank])^2 ## need to take square to get eigenvalues
      Um <- rpm$u[,1:rank]
      ## correct the sign (such that the eigenvectors are orthonormal)
      signflip = which(diag(t(Um)%*%Um)<0)
      Um[,signflip] = -Um[,signflip]
    }
    invSIG <- (diag(1,n)-Um%*%diag(1/((1/(tau2_nonadditive*Dm))+1))%*%t(Um))/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(tau2_nonadditive*Dm+1)) ## includes sigma2

  }else if(invmethod=="lzSMWmax"){ ## sherman morrison woodbury

    if(only_tau2==TRUE & !is.null(Um) & !is.null(Dm)){ ## if only updating tau2, we can skip the approximation and just manipulate tau2
      Um <- Um
      Dm <- Dm
    }else{
      lz <- mgcv::slanczos(eigenQuadProd(K$K,P),k=rank)#mgcv::slanczos(t(P)%*%K$K%*%P,k=rank) ## highest m eigenvalues
      Um <- lz$vectors
      Dm <- lz$values # does not include tau2
    }
    invSIG <- (diag(1,n)-Um%*%diag(1/((1/(tau2_nonadditive*Dm))+1))%*%t(Um))/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(tau2_nonadditive*Dm+1)) ## includes sigma2

  }else if(invmethod=="lzSMWmin"){

    if(only_tau2==TRUE & !is.null(Um) & !is.null(Dm)){ ## if only updating tau2, we can skip the approximation and just manipulate tau2
      Um <- Um
      Dm <- Dm
    }else{
      lz <- mgcv::slanczos(eigenQuadProd(K$K,P),k=0,kl=rank) #mgcv::slanczos(t(P)%*%K$K%*%P,k=0,kl=rank) ## lowest m eigenvalues
      Um <- lz$vectors
      Dm <- lz$values # does not include tau2
    }
    invSIG <- (diag(1,n)-Um%*%diag(1/((1/(tau2_nonadditive*Dm))+1))%*%t(Um))/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(tau2_nonadditive*Dm+1)) ## includes sigma2

  }else if(invmethod=="lzSMWboth"){

    if(only_tau2==TRUE & !is.null(Um) & !is.null(Dm)){ ## if only updating tau2, we can skip the approximation and just manipulate tau2
      Um <- Um
      Dm <- Dm
    }else{
      lz <- mgcv::slanczos(eigenQuadProd(K$K,P),k=rank,kl=rank) #mgcv::slanczos(t(P)%*%K$K%*%P,k=rank,kl=rank) ## highest m and lowest m eigenvalues
      Um <- lz$vectors
      Dm <- lz$values # does not include tau2
    }
    invSIG <- (diag(1,n)-Um%*%diag(1/((1/(tau2_nonadditive*Dm))+1))%*%t(Um))/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(tau2_nonadditive*Dm+1)) ## includes sigma2

  }else if(invmethod=="lzDIRECTmax"){

    lz <- mgcv::slanczos(diag(1,n)+tau2_nonadditive*eigenQuadProd(K$K,P),k=rank) #mgcv::slanczos(diag(1,n)+tau2_nonadditive*t(P)%*%K$K%*%P,k=rank) ## highest m eigenvalues
    Um <- lz$vectors
    Dm <- lz$values # does not include tau2
    invSIG <- Um%*%diag(1/Dm)%*%t(Um)/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(Dm)) ## includes sigma2

  }else if(invmethod=="lzDIRECTmin"){

    lz <- mgcv::slanczos(diag(1,n)+tau2_nonadditive*eigenQuadProd(K$K,P),k=rank,kl=rank) #mgcv::slanczos(diag(1,n)+tau2_nonadditive*t(P)%*%K$K%*%P,k=rank,kl=rank) ## lowest m eigenvalues
    Um <- lz$vectors
    Dm <- lz$values # does not include tau2
    invSIG <- Um%*%diag(1/Dm)%*%t(Um)/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(Dm)) ## includes sigma2

  }else if(invmethod=="lzDIRECTboth"){

    lz <- mgcv::slanczos(diag(1,n)+tau2_nonadditive*eigenQuadProd(K$K,P),k=rank,kl=rank) #mgcv::slanczos(diag(1,n)+tau2_nonadditive*t(P)%*%K$K%*%P,k=rank,kl=rank) ## highest m and lowest m eigenvalues
    Um <- lz$vectors
    Dm <- lz$values # does not include tau2
    invSIG <- Um%*%diag(1/Dm)%*%t(Um)/sigma2 ##SMW formula ## includes sigma2
    logdetSIG <- n*log(sigma2)+sum(log(Dm)) ## includes sigma2

  }

  return(list(invSIG=invSIG,
              logdetSIG=logdetSIG,
              Um=Um,
              Dm=Dm))
}
