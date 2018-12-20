#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <wishart.h>
#include <mvnorm.h>
#include <omp.h>

using namespace arma;
using namespace std;
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
bool kMeans(mat &C, cube &P, unsigned int K,int n_iter){
  uword d=P.n_slices;
  uword N=P.n_cols*P.n_rows;
  mat data(d,N);
  vec temp;
  for(unsigned int i=0;i<P.n_rows;i++)
    for(unsigned int j=0;j<P.n_cols;j++){
      temp=P.tube(i,j);
      data.col(i*P.n_rows+j)=temp;
    }
    
    mat means;
    bool status = kmeans(means, data, K, random_subset, n_iter, false);
    int Minx=1000000;
    vec temp1,temp2;
    for(unsigned int i=0;i<P.n_rows;i++)
      for(unsigned int j=0;j<P.n_cols;j++){
        temp1=P.tube(i,j);
        Minx=1000000;
        for(unsigned int k=0;k<K;k++){
          temp2=means.col(k);
          if(norm(temp1-temp2)<Minx)
            C(i,j)=k+1;
        }
      }
        
 return status; 
}



// [[Rcpp::export]]
mat ImageGibbs(int K, cube P, mat C, mat MuK, cube SigmaK, double alpha, double beta, 
                 vec mu0, mat lambda0, vec v0, mat sigma0, int burnIn,int mcmcN){
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
  unsigned int i,j,k;
  int m=P.n_rows,n=P.n_cols;
  cube sampleC(m,n,mcmcN);
  cube Sigma(3,3,K);
  mat Mu(3,K);
  Sigma=SigmaK;
  Mu=MuK;
  
  mat mu1(3,1);
  mat sigma1(3,3);
  //Update
  for(i=0;i<burnIn+ mcmcN;i++){
    
    // update Mu
    for(k=0;k<K;k++){
      
    }
    
    
    
    
    
    
  }
  
  return sampleC;
  
}

/*** R
# read in image
library(jpeg)
P = readJPEG("./dog.jpg",FALSE)
C = P[,,1]
K=15
n_iter=10000;
kMeans(C,P,K,n_iter)
*/
