// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//


// [[Rcpp::export]]
arma::mat solvecpp(const arma::mat& A){
    return(pinv(A, 1.0e-25));}




// [[Rcpp::export]]
double maxeigencpp( const arma::mat& X) {
    return(max(eig_sym( X )));
}


// [[Rcpp::export]]
double normcpp(const arma::vec& x, const arma::mat& A, bool euc=true) {
    if (euc == true)  {return(norm(x,2));}
    else {return(as_scalar(sqrt(trans(x) *A* x) ));}
}
