#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>

using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//Potts Model
//compute cij
inline int mirrorIndex(int fetchI, int length)
{
  return (fetchI < 0) + (fetchI >= length) * (length - 2);
}
// [[Rcpp::export]]
double PrCij(mat &C, double alpha, double beta, int i, int j)
{
  uword m = C.n_rows, n = C.n_cols;
  int cij = C.at(i, j);
  double sum = 0;
  sum += (C.at(mirrorIndex(i - 1, m), j) == cij ? alpha : beta);
  sum += (C.at(mirrorIndex(i + 1, m), j) == cij ? alpha : beta);
  sum += (C.at(i, mirrorIndex(j + 1, n)) == cij ? alpha : beta);
  sum += (C.at(i, mirrorIndex(j - 1, n)) == cij ? alpha : beta);
  return sum;
}

// [[Rcpp::export]]
double Pr(mat &C, double alpha, double beta)
{
  double sum = 0;
  for (unsigned int i = 0; i < C.n_rows; i++)
  {
    for (unsigned int j = 0; j < C.n_cols; j++)
    {
      sum += PrCij(C, alpha, beta, i, j);
    }
  }
  return sum;
}

// [[Rcpp::export]]
double Pr_parallel(mat &C, double alpha, double beta)
{
  double sum = 0;
  uword i, j;
#pragma omp parallel for reduction(+                                                            \
                                   : sum) schedule(static) private(i, j) shared(C, alpha, beta) \
    collapse(2)
  for (i = 0; i < C.n_rows; i++)
  {
    for (j = 0; j < C.n_cols; j++)
    {
      sum += PrCij(C, alpha, beta, i, j);
    }
  }
  return sum;
}
/*** R
library(microbenchmark)
microbenchmark(Pr_parallel(C, alpha, beta),Pr(C, alpha, beta))
*/