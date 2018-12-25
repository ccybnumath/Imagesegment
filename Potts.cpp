#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#include <mvnorm.h>
#include <omp.h>

using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//Potts Model
//compute cij
int mirrorIndex(int fetchI, int length){
  return (fetchI<0)+(fetchI>=length)*(length-2);
}

double PrCij(mat &C, double alpha, double beta, int i, int j, int cij){
  vector<int> a{-1,1};
  unsigned int m=0;
  double sum=0;
  for(m=0;m<1;m++){
    sum+=(C.at(mirrorIndex(i+a.at(m),C.n_rows),j)==cij?alpha:beta);
    sum+=(C.at(i, mirrorIndex(j+a.at(m),C.n_cols))==cij?alpha:beta);
  }
  return sum;
}

// [[Rcpp::export]]
double Pr(mat &C, double alpha, double beta){
  double sum=0;
  for(unsigned int i=0;i<C.n_rows;i++){
    for(unsigned int j=0;j<C.n_cols;j++){
      sum+=PrCij(C,alpha,beta,i,j,C.at(i,j));
    }
  }
  return sum;
}