#include <RcppArmadillo.h>
using namespace Rcpp;


// do_rgiglocal uses rgig which is imported from GIGrvg
// [[Rcpp::depends(GIGrvg,RcppArmadillo)]]
double do_rgiglocal(double lambda, double chi, double psi) {
  
  SEXP (*fun)(int,double,double,double) = NULL;
  
  if (!fun)
    fun = (SEXP(*)(int,double,double,double)) R_GetCCallable("GIGrvg", "do_rgig");
  
  return as<double>(fun(1,lambda,chi,psi));
}