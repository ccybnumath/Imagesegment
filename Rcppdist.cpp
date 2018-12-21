#include <RcppDist.h> 
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace arma;
using namespace std;

// [[Rcpp::export]]
mat iwish(int df, mat S){
  return riwish(df,S);
}

// [[Rcpp::export]]
vec f(){
  cube M(3,4,5,fill::randn);
  vec temp;
  temp=M.tube(1,1);
  return temp;
}

vec sumP(cube &P, mat &C, int k){
  vec temp;
  vec sum(3,fill::zeros);
  for(unsigned i=0;i<C.n_rows;i++)
    for(unsigned j=0;j<C.n_cols;j++){
      if(C.at(i,j)==k){
        temp=P.tube(i,j);
        sum+=temp;
      }
    }
    return sum;
}


// [[Rcpp::export]]
vec ImageGibbs(unsigned int K, cube P, mat C, mat Mu, cube Sigma, double alpha, double beta, 
               vec mu0, mat lambda0, int v0, mat sigma0, unsigned int burnIn,unsigned int mcmcN){
  vec temp;
  temp=sumP(P,C,K);
  return temp;
}
   