#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List rp(int n, int k, int alpha, arma::mat C , int m){
  // in the C code, we have to define character of value

  mat A(n, k, fill::randn);
  mat B = A * pow(k, -0.25) ;
  mat phi = pow(C, alpha) * B;
  mat aK = phi.t() * C * phi;
  //
  mat U;
  vec s;
  mat V;
  svd(U, s, V, aK) ;
  mat H= diagmat(pow(s,-0.5));
  mat Km= (C * phi)*U*(H);
  //Eigen Decomposition
  mat Ub;
  vec sb;
  mat Vb;
  svd(Ub, sb, Vb, Km);

  return List::create(Named("u") = Ub, Named("d") = sb );
}


// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::export]]
SEXP eigenQuadProd(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = B.transpose() * A * B;

  return Rcpp::wrap(C);
}
