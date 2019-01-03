#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <mvnorm.h>
#include <omp.h>
#include <wishart.h>

#include <UpdateC.cpp>
using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
//User-defined reduction
#pragma omp declare reduction(+                    \
                              : arma::mat          \
                              : omp_out += omp_in) \
    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(+                    \
                              : arma::vec          \
                              : omp_out += omp_in) \
    initializer(omp_priv = omp_orig)

// [[Rcpp::export]]
mat computeSk_parallel(cube &P, mat &C, mat &Mu, int k)
{
  mat Sk(3, 3, fill::zeros);
  vec temp;
  uword i, j;
#pragma omp parallel for schedule(static) reduction(+                                             \
                                                    : Sk) private(i, j, temp) shared(C, P, Mu, k) \
    collapse(2)
  for (i = 0; i < C.n_rows; i++)
    for (j = 0; j < C.n_cols; j++)
    {
      if (C.at(i, j) == k)
      {
        temp = P.tube(i, j);
        temp -= Mu.col(k);
        Sk += temp * temp.t();
      }
    }

  return Sk;
}

// [[Rcpp::export]]
mat computeSk(cube &P, mat &C, mat &Mu, int k)
{
  mat Sk(3, 3, fill::zeros);
  vec temp;
  uword i, j;
  for (i = 0; i < C.n_rows; i++)
    for (j = 0; j < C.n_cols; j++)
    {
      if (C.at(i, j) == k)
      {
        temp = P.tube(i, j);
        temp -= Mu.col(k);
        Sk += temp * temp.t();
      }
    }

  return Sk;
}

// [[Rcpp::export]]
vec sumP_parallel(cube &P, mat &C, int k)
{
  vec temp;
  vec sum1(3, fill::zeros);
  uword i, j;
#pragma omp parallel for schedule(static) reduction(+                                           \
                                                    : sum1) private(i, j, temp) shared(C, P, k) \
    collapse(2)
  for (i = 0; i < C.n_rows; i++)
    for (j = 0; j < C.n_cols; j++)
    {
      if (C.at(i, j) == k)
      {
        temp = P.tube(i, j);
        sum1 += temp;
      }
    }
  return sum1;
}

// [[Rcpp::export]]
vec sumP(cube &P, mat &C, int k)
{
  vec temp;
  vec sum1(3, fill::zeros);
  for (unsigned i = 0; i < C.n_rows; i++)
    for (unsigned j = 0; j < C.n_cols; j++)
    {
      if (C.at(i, j) == k)
      {
        temp = P.tube(i, j);
        sum1 += temp;
      }
    }
  return sum1;
}

// [[Rcpp::export]]
cube ImageGibbs(uword K, cube P, mat &C, mat &Mu, cube &Sigma, double alpha, double beta,
                vec mu0, mat lambda0, int v0, mat sigma0, uword burnIn, uword mcmcN)
{
  /*
   * * INPUT:
   * K represents K groups
   * P is the RGB image Cube
   * m,n are the rows & cols of the Image P
   * C represents the membership of the Image,cij belongs to {1,2,...,K}
   * alpha, beta are the parameters of Potts Models
   * mu0 is the mean of the prior Normal distribution 
   * lambda0 is the variance Matrix
   * v0, sigma0 is the parameters of the Inverse-Wishart
   * * OUTPUT
   * sample
   */
  //Init
  uword k, l;
  uword m = P.n_rows, n = P.n_cols;
  cube sampleC(m, n, mcmcN);
  C = C - 1;

  mat tempMu(3, 1);
  mat tempSigma(3, 3);
  int tempnk = 0;

  int tempv0 = 0;
  mat tempS(3, 3, fill::zeros);
  mat tempSk(3, 3, fill::zeros);
  //Update
  for (l = 0; l < burnIn + mcmcN; l++)
  {

    // update Mu
    for (k = 0; k < K; k++)
    {
      tempnk = sum(sum(C == k));
      tempSigma = inv(tempnk * inv(Sigma.slice(k)) + inv(lambda0));
      tempMu = tempSigma * (inv(Sigma.slice(k)) * sumP_parallel(P, C, k) + inv(lambda0) * mu0);
      Mu.col(k) = rmvnorm(1, tempMu, tempSigma).t();
    }

    //update Sigma
    for (k = 0; k < K; k++)
    {
      tempnk = sum(sum(C == k));
      tempv0 = v0 + tempnk;
      tempSk = computeSk_parallel(P, C, Mu, k);
      tempSk = inv(sigma0) + tempSk;
      Sigma.slice(k) = riwish(tempv0, tempSk);
    }

    //update Cij
    mat Prob(m, n, fill::zeros);
    InitProb_parallel(Prob, C, alpha, beta);
    UpdateC_parallel(C, P, Mu, Sigma, m, n, K, alpha, beta, Prob);

    //record
    if (l >= burnIn)
      sampleC.slice(l - burnIn) = C;
  }

  return sampleC;
  // to do find Mod of sample C;
}
