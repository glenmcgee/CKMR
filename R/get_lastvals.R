
#' Get last values.
#'
#' Extract parameter values from last iteration.
#'
#' @param fit Fitted model
#'
#' @return Returns list of last values
#'
#' @export
get_lastvals <- function(fit){
  n_samples <- nrow(fit$theta)
  return(list(theta=fit$theta[n_samples,],
              pi_additive=fit$pi_additive[n_samples],
              gamma_additive=fit$gamma_additive[n_samples,],
              beta_additive=fit$beta_additive[n_samples,],
              tau2_additive=fit$tau2_additive[n_samples,],
              pi_nonadditive=fit$pi_nonadditive[n_samples],
              gamma_nonadditive=fit$gamma_nonadditive[n_samples,],
              rho_nonadditive =fit$rho_nonadditive[n_samples,],
              tau2_nonadditive=fit$tau2_nonadditive[n_samples],
              alpha=fit$alpha[n_samples,],
              sigma2=fit$sigma2[n_samples],
              SS=fit$SS,
              B=fit$B,
              P=fit$P))
}




