#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void print(){
  cout<<"Hello from Rcpp.";
}