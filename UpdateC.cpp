#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#include <mvnorm.h>
#include <omp.h>

#include <Potts.cpp>
using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//Update Cij
// [[Rcpp::export]]
int UpdateCij(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, 
              uword n, uword K, double alpha, double beta, 
              uword i, uword j){
  uword k;
  vec probK(K,fill::zeros);
  vec N(K,fill::zeros);
  vec fullvec = regspace<vec>(0,K-1);
  mat Cij=C;
  for(k=0;k<K;k++){
    Cij.at(i,j)=k;
    N.at(k)=sum(dmvnorm(P.tube(i,j),Mu.col(k),Sigma.slice(k)));
    probK.at(k)=Pr(Cij,alpha,beta);
  }
  probK-=max(probK);
  probK=exp(probK);
  N/=max(N);
  probK=probK%N/sum(probK.t()*N);
  return sum(Rcpp::RcppArmadillo::sample(fullvec,1,true,probK));
}

// [[Rcpp::export]]
int UpdateCij_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, 
              uword n, uword K, double alpha, double beta, 
              uword i, uword j){
  uword k;
  vec probK(K,fill::zeros);
  vec N(K,fill::zeros);
  vec fullvec = regspace<vec>(0,K-1);
  mat Cij=C;
  for(k=0;k<K;k++){
    Cij.at(i,j)=k;
    N.at(k)=sum(dmvnorm(P.tube(i,j),Mu.col(k),Sigma.slice(k)));
    probK.at(k)=Pr_parallel(Cij,alpha,beta);
  }
  probK-=max(probK);
  probK=exp(probK);
  N/=max(N);
  probK=probK%N/sum(probK.t()*N);
  return sum(Rcpp::RcppArmadillo::sample(fullvec,1,true,probK));
}

// [[Rcpp::export]]
void UpdateC(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta){
  uword i,j;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++){
      C.at(i,j)=UpdateCij(C,P,Mu,Sigma,m,n,K,alpha,beta,i,j);
    }
    
}

// [[Rcpp::export]]
void UpdateC_parallel(mat &C, cube &P, mat &Mu, cube &Sigma, uword m, uword n, uword K, double alpha, double beta){
  uword i,j;
#pragma omp parallel for schedule(static)            \
  private(i,j) shared(m,n,C,P,Mu,Sigma,K,alpha,beta) \
    collapse(2)
    for(i=0;i<m;i++)
      for(j=0;j<n;j++){
        C.at(i,j)=UpdateCij(C,P,Mu,Sigma,m,n,K,alpha,beta,i,j);
      }
      
}
/*** R
library(microbenchmark)
microbenchmark(UpdateCij(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19),UpdateCij_parallel(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19))
*/