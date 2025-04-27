
#' Get B and P matrices
#'
#' @keywords internal
get_BP <- function(xtheta,B,SS,j,d,n,approxProj,P=NULL){

  ## update only the relevant columns of B (corresponding the the jth index)
  B[,sum(d[0:(j-1)])+(1:d[j])] <- mgcv::PredictMat(SS[[j]],data=xtheta)
  Bproj <- cbind(1,B[,apply(B,2, function(x) stats::var(x)!=0)]) ## remove columns that are constant due to gamma_additive=0 (hence non invertible)
  # P <- diag(1,n)-Bproj%*%solve(t(Bproj)%*%Bproj,t(Bproj))
  if(approxProj==FALSE | is.null(P)){ ## compute new P, but only if we are not using the approximate version
    # P <- diag(1,n)-eigenMapMatMult(Bproj, solve(eigenMapMatMult(t(Bproj),Bproj),t(Bproj))) ## using eigen. could use for solve as well...
    try(invBB <- solve(eigenMapMatMult(t(Bproj),Bproj)),
    invBB <- MASS::ginv(eigenMapMatMult(t(Bproj),Bproj)))
    P <- diag(1,n)-eigenMapMatMult(Bproj,eigenMapMatMult(invBB,t(Bproj))) ##using generalized inverse due to multicollinearity
  }

  return(list(B=B,
              P=P))
}
