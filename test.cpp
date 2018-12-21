#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#include <mvnorm.h>
#include <omp.h>

using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

int mirrorIndex(int fetchI, int length){
  return (fetchI<0)+(fetchI>=length)*(length-2);
}


//Potts Model
double Pr(mat &C, double alpha, double beta, int i, int j, int k){
  vector<int> a{-1,1};
  int m=0;
  double sum=0;
  for(m=0;m<1;m++)
    sum+=(C.at(mirrorIndex(i+a.at(m),C.n_rows),j)==k?alpha:beta);
  return exp(sum);
}


mat computeSk(cube &P, mat &C, mat &Mu, int k){
  mat Sk(3,3,fill::zeros);
  vec temp;
  for(unsigned i=0;i<C.n_rows;i++)
    for(unsigned j=0;j<C.n_cols;j++){
      if(C.at(i,j)==k){
        temp=P.tube(i,j);
        temp-=Mu.col(k);
        Sk+=temp*temp.t();
      }
    }
    
    return Sk;
}

vec sumP(cube &P, mat &C, int k){
  vec temp;
  vec sum1(3,fill::zeros);
  for(unsigned i=0;i<C.n_rows;i++)
    for(unsigned j=0;j<C.n_cols;j++){
      if(C.at(i,j)==k){
        temp=P.tube(i,j);
        sum1+=temp;
      }
    }
    return sum1;
}

// [[Rcpp::export]]
cube ImageGibbs(unsigned int K, cube P, mat C, mat Mu, cube Sigma, double alpha, double beta, 
                vec mu0, mat lambda0, int v0, mat sigma0, unsigned int burnIn,unsigned int mcmcN){
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
  unsigned int i,j,k,l;
  unsigned int m=P.n_rows,n=P.n_cols;
  cube sampleC(m,n,mcmcN);
  C=C-1;
  
  mat tempMu(3,1);
  mat tempSigma(3,3);
  int tempnk=0;
  
  int tempv0=0;
  mat tempS(3,3,fill::zeros);
  mat tempSk(3,3,fill::zeros);
  //Update
  for(l=0;l<burnIn+ mcmcN;l++){
    
    // update Mu
    for(k=0;k<K;k++){
      tempnk=sum(sum(C==k));
      tempSigma=inv(tempnk*inv(Sigma.slice(k))+inv(lambda0));
      tempMu=tempSigma*(inv(Sigma.slice(k))*sumP(P,C,k)+inv(lambda0)*mu0);
      Mu.col(k)=rmvnorm(1,tempMu,tempSigma).t();
    }
    
  }
  
  return sampleC;
  
}
