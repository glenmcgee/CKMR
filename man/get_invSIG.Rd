% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_invSIG.R
\name{get_invSIG}
\alias{get_invSIG}
\title{Get Sigma inverse and log of det(Sigma)}
\usage{
get_invSIG(
  P,
  K,
  n = nrow(P),
  tau2_nonadditive,
  sigma2,
  invmethod = "exact",
  rank = NULL,
  only_tau2 = FALSE,
  Um = NULL,
  Dm = NULL
)
}
\description{
Get Sigma inverse and log of det(Sigma)
}
\keyword{internal}
