#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>

using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//Update Cij
// to do parallel
// [[Rcpp::export]]
void dparallel()
{
  int i, j;
  int m = 3, n = 3;
  int total = omp_get_max_threads();
#pragma omp parallel for schedule(static) private(i) \
    collapse(2)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (j % 2 == (i + 1) % 2)
      {
        printf("Location: %d,%d\n", i, j);
      }
      //printf("Hello from location %d, %d\n",i,j);
    }
  }
}
/*** R
dparallel()
  */