# require(refund)

#' Get basis
#'
#' @keywords internal
getbasis <- function(X, ## design matrix for index (columns represent different time points)
                     basis.opts){ ## choice of basis


  #account for missing input in basis function
  if(toupper(basis.opts$type)=="NONE"){

    return(list(psi=diag(ncol(X)), type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="RANKED"){

    return(list(psi=upper.tri(diag(ncol(X)),diag=TRUE)+0, type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="PCA"){

    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) basis.opts$pve <- .99
    }
    pca <- stats::prcomp(X)
    nbases <- min(which(summary(pca)$importance[3,]>basis.opts$pve)) ## number of PCs to include
    return(list(psi=as.matrix(pca$rotation[,1:nbases]), eigenvalues=pca$sdev[1:nbases]^2, pve=basis.opts$pve,type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="NS"){

    if(is.null(basis.opts$df) & !is.null(basis.opts$knots)){
      basis.opts$df <- basis.opts$knots+2
    }else if(is.null(basis.opts$df)){
      basis.opts$df <- 5
    }

    B1 <- splines::ns(seq(1,ncol(X)),df=basis.opts$df, intercept=TRUE)
    X <-  B1 %*% qr.solve(B1,t(X))
    svdX <- svd(X)

    return(list(psi=svdX$u[,1:ncol(B1)], lambda=svdX$d[1:ncol(B1)], pve=basis.opts$pve,type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="BS"){

    if(is.null(basis.opts$df)) basis.opts$df <- round(ncol(X)/5)
    B1 <- splines::bs(seq(1,ncol(X)),df=basis.opts$df, intercept=TRUE)
    X <-  B1 %*% qr.solve(B1,t(X))
    svdX <- svd(X)

    return(list(psi=svdX$u[,1:ncol(B1)], lambda=svdX$d[1:ncol(B1)], pve=basis.opts$pve,type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="FACE"){

    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) pve <- .99
    }

    if(is.null(basis.opts$knots)) basis.opts$knots <- round(ncol(X)/3)
    fitted <- refund::fpca.face(X, knots=basis.opts$knots,pve=basis.opts$pve)
    lambda <- fitted$evalues
    return(list(psi=fitted$efunctions, lambda=lambda, knots=basis.opts$knots,pve=basis.opts$pve,type=basis.opts$type))

  }else if(toupper(basis.opts$type)=="GAM"){

    #correct pve if needed.
    if(is.null(basis.opts$pve) | !is.numeric(basis.opts$pve)){
      basis.opts$pve=.99
    }else{
      if(basis.opts$pve<0 | basis.opts$pve>1) pve <- .99
    }


    for(i in 1:nrow(X)){
      X[i,] <- mgcv::predict.gam(mgcv::gam(X[i,]~mgcv::s(seq(1:ncol(X)))))
    }

    eigcorX <- eigen(stats::cor(X))
    nbases <- min(which(cumsum(eigcorX$values)/sum(eigcorX$values)>basis.opts$pve))
    return(list(psi=eigcorX$vectors[,1:nbases], eigenvalues=eigcorX$values[1:nbases], pve=basis.opts$pve, df=basis.opts$df,type=basis.opts$type))


  }else{

    stop("basis type not recognized.")

  }


}
