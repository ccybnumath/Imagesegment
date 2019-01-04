#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>



using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]

// [[Rcpp::export]]
void f(cube &P,mat &Mu,cube &Sigma){
  int i,j;
  int m = P.n_rows;
  int n = P.n_cols;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n,P,Mu,Sigma)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      dmvnorm(P.tube(i, j), Mu.col(1), Sigma.slice(1));
    }
  }
}

// [[Rcpp::export]]
void g(cube &P,mat &Mu,cube &Sigma){
 double x;
  int i,j;
  int m = P.n_rows;
  int n = P.n_cols;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n,P,Mu,Sigma)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      x=sum(dmvnorm(P.tube(i, j), Mu.col(1), Sigma.slice(1)));
    }
  }
}
/*** R
library(microbenchmark)
microbenchmark(f(P,Mu,Sigma),g(P,Mu,Sigma),times = 1)
*/

