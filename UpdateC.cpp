#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <trng/discrete_dist.hpp>
#include <trng/yarn2.hpp>
#include <wishart.h>

#include <Potts.cpp>
using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]

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
int UpdateCij_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m,
                       uword n, uword K, double alpha, double beta,
                       uword i, uword j, mat Prob)
{
  uword k;
  vec probK(K, fill::zeros);
  vec N(K, fill::zeros);
  vec fullvec = regspace<vec>(0, K - 1);
  mat Cij = C;
  mat Probij = Prob;

  //trng distSampling
  trng::yarn2 rx;
  double x;
  rx.seed(10);
  int size = omp_get_num_threads(); // get total number of processes
  int rank = omp_get_thread_num();  // get rank of current process
  rx.split(size, rank);             // choose sub-stream no. rank out of size streams

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
  //to do
  //Normalization
  probK -= max(probK);
  probK = exp(probK);
  probK /= max(probK);
  N /= max(N);
  //cout << N << endl;
  probK = probK % N / sum(probK.t() * N);
  //cout << probK << endl;
  trng::discrete_dist distSampling(probK.begin(), probK.end());
  return distSampling(rx);
  // return sum(Rcpp::RcppArmadillo::sample(fullvec, 1, true, probK));
}

// [[Rcpp::export]]
void UpdateC_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta, mat Prob)
{

  uword i, j;
#pragma omp parallel for schedule(static) private(i, j) shared(m, n, C, P, Mu, Sigma, K, alpha, beta, Prob) \
    collapse(2)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      if (i % 2 == j % 2)
        C.at(i, j) = UpdateCij_parallel(C, P, Mu, Sigma, m, n, K, alpha, beta, i, j, Prob);
    }
  }
#pragma omp parallel for schedule(static) private(i, j) shared(m, n, C, P, Mu, Sigma, K, alpha, beta, Prob) \
    collapse(2)
  for (i = 0; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      if ((i + 1) % 2 == j % 2)
        C.at(i, j) = UpdateCij_parallel(C, P, Mu, Sigma, m, n, K, alpha, beta, i, j, Prob);
    }
  }
}

/*
  * 我的电脑的结果
  * > microbenchmark(UpdateC_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
  *Unit: seconds
  *                                                                   expr      min       lq     mean   median       uq      max neval
  * UpdateC_parallelUnique(C, P, Mu, Sigma, m, n, K, alpha, beta,      Prob) 451.9796 451.9796 451.9796 451.9796 451.9796 451.9796     1
  * 
  * 服务器上面的
  * microbenchmark(UpdateC_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
  * Unit: seconds
  *                                                                     expr      min       lq     mean   median       uq      max neval
  *  UpdateC_parallelUnique(C, P, Mu, Sigma, m, n, K, alpha, beta,      Prob) 22.18771 22.18771 22.18771 22.18771 22.18771 22.18771     1
  */