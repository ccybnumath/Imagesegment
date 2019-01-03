#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>

#include <UpdateC.cpp>
using namespace arma;
using namespace std;

// [[Rcpp::export]]

