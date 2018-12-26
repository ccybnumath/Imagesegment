#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>

#include <Potts.cpp>
using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//Update Cij
// to do parallel
// [[Rcpp::export]]
int UpdateCij(mat &C, cube &P, mat &Mu, cube &Sigma, uword m,
              uword n, uword K, double alpha, double beta,
              uword i, uword j)
{
  uword k;
  vec probK(K, fill::zeros);
  vec N(K, fill::zeros);
  vec fullvec = regspace<vec>(0, K - 1);
  mat Cij = C;
  for (k = 0; k < K; k++)
  {
    Cij.at(i, j) = k;
    N.at(k) = sum(dmvnorm(P.tube(i, j), Mu.col(k), Sigma.slice(k)));
    probK.at(k) = Pr(Cij, alpha, beta);
  }
  probK -= max(probK);
  probK = exp(probK);
  probK /= max(probK);
  //cout << probK << endl;
  N /= max(N);
  //cout << N << endl;
  probK = probK % N / sum(probK.t() * N);
  return sum(Rcpp::RcppArmadillo::sample(fullvec, 1, true, probK));
}

// [[Rcpp::export]]
int UpdateCij_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m,
                       uword n, uword K, double alpha, double beta,
                       uword i, uword j)
{
  uword k;
  vec probK(K, fill::zeros);
  vec N(K, fill::zeros);
  vec fullvec = regspace<vec>(0, K - 1);
  mat Cij = C;
  for (k = 0; k < K; k++)
  {
    Cij.at(i, j) = k;
    N.at(k) = sum(dmvnorm(P.tube(i, j), Mu.col(k), Sigma.slice(k)));
    probK.at(k) = Pr_parallel(Cij, alpha, beta);
  }
  probK -= max(probK);
  probK = exp(probK);
  probK /= max(probK);
  //cout << probK << endl;
  N /= max(N);
  //cout << N << endl;
  probK = probK % N / sum(probK.t() * N);
  return sum(Rcpp::RcppArmadillo::sample(fullvec, 1, true, probK));
}

// [[Rcpp::export]]
void UpdateC(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta)
{
  uword i, j;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      C.at(i, j) = UpdateCij(C, P, Mu, Sigma, m, n, K, alpha, beta, i, j);
    }
}

// [[Rcpp::export]]
void UpdateC_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta)
{
  uword i, j;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n, C, P, Mu, Sigma, K, alpha, beta) \
    collapse(2)
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      C.at(i, j) = UpdateCij(C, P, Mu, Sigma, m, n, K, alpha, beta, i, j);
    }
}

// [[Rcpp::export]]
void InitProb_parallel(mat &Prob, mat &C, double alpha, double beta)
{
  //Init Prob
  uword i, j, m, n;
  m = C.n_rows;
  n = C.n_cols;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n, C, alpha, beta, Prob)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      Prob.at(i, j) = PrCij(C, alpha, beta, i, j);
    }
  }
}

// [[Rcpp::export]]
int UpdateCij_parallelUnique(mat &C, cube &P, mat &Mu, cube &Sigma, uword m,
                             uword n, uword K, double alpha, double beta,
                             uword i, uword j, mat Prob)
{
  uword k;
  vec probK(K, fill::zeros);
  vec N(K, fill::zeros);
  vec fullvec = regspace<vec>(0, K - 1);
  mat Cij = C;
  mat Probij = Prob;
  for (k = 0; k < K; k++)
  {
    //revise Probij
    Cij.at(i, j) = k;
    Probij.at(i, j) = PrCij(Cij, alpha, beta, i, j);
    Probij.at(mirrorIndex(i - 1, m), j) = PrCij(Cij, alpha, beta, mirrorIndex(i - 1, m), j);
    Probij.at(mirrorIndex(i + 1, m), j) = PrCij(Cij, alpha, beta, mirrorIndex(i + 1, m), j);
    Probij.at(i, mirrorIndex(j - 1, n)) = PrCij(Cij, alpha, beta, i, mirrorIndex(j - 1, n));
    Probij.at(i, mirrorIndex(j + 1, n)) = PrCij(Cij, alpha, beta, i, mirrorIndex(j + 1, n));

    N.at(k) = sum(dmvnorm(P.tube(i, j), Mu.col(k), Sigma.slice(k)));
    probK.at(k) = accu(Probij);
  }

  //Normalization
  probK -= max(probK);
  probK = exp(probK);
  probK /= max(probK);
  N /= max(N);
  //cout << N << endl;
  probK = probK % N / sum(probK.t() * N);
  //cout << probK << endl;
  return sum(Rcpp::RcppArmadillo::sample(fullvec, 1, true, probK));
}

// [[Rcpp::export]]
void UpdateC_parallelUnique(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta, mat Prob)
{
  /*
  * 我的电脑的结果
  * > microbenchmark(UpdateC_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
  *Unit: seconds
  *                                                                   expr      min       lq     mean   median       uq      max neval
  * UpdateC_parallelUnique(C, P, Mu, Sigma, m, n, K, alpha, beta,      Prob) 451.9796 451.9796 451.9796 451.9796 451.9796 451.9796     1
  * 
  * 服务器上面的
  * > microbenchmark(UpdateC_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
  *Unit: seconds
  *                                                                   expr      min       lq     mean   median       uq
  * UpdateC_parallelUnique(C, P, Mu, Sigma, m, n, K, alpha, beta,      Prob) 25.82491 25.82491 25.82491 25.82491 25.82491
  *     max neval
  * 25.82491     1 */
  */
  uword i, j;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n, C, P, Mu, Sigma, K, alpha, beta, Prob) \
    collapse(2)
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
    {
      C.at(i, j) = UpdateCij_parallelUnique(C, P, Mu, Sigma, m, n, K, alpha, beta, i, j, Prob);
    }

}

/*** R
library(microbenchmark)
microbenchmark(UpdateCij(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19),UpdateCij_parallel(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19),UpdateCij_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19,Prob))
microbenchmark(UpdateC_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
*/