#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#include <mvnorm.h>
#include <omp.h>

using namespace arma;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
int f(int K){
  uword d = 5;       // dimensionality
  uword N = 10000;   // number of vectors
  
  mat data(d, N, fill::randu);
  
  mat means;
  
  bool status = kmeans(means, data, K, random_subset, 10, true);
  
  if(status == false)
  {
    cout << "clustering failed" << endl;
  }
  
  means.print("means:");
  return 1;
}