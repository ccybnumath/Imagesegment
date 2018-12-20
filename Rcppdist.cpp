#include <RcppDist.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat iwish(int df, mat S){
  return riwish(df,S);
}

// [[Rcpp::export]]
vec f(){
  cube M(3,4,5,fill::randn);
  vec temp;
  temp=M.tube(1,1);
  return temp;
}
