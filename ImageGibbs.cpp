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
      if(C.at(i,j)==k)
        temp=P.tube(i,j);
        temp-=Mu.col(k);
        Sk+=temp*temp.t();
    }
  
  return Sk;
}



// [[Rcpp::export]]
mat ImageGibbs(int K, cube P, mat C, mat Mu, cube Sigma, double alpha, double beta, 
                 vec mu0, mat lambda0, int v0, mat sigma0, int burnIn,int mcmcN){
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
  int m=P.n_rows,n=P.n_cols;
  cube sampleC(m,n,mcmcN);
  C=C-1;
  
  mat tempMu(3,1);
  mat tempSigma(3,3);
  int tempnk=0;
  
  int tempv0=0;
  mat tempS(3,3,fill::zeros);
  mat tempSk(3,3,fill::zeros);
  //Update
  for(l=0;i<burnIn+ mcmcN;i++){
    
    // update Mu
    for(k=0;k<K;k++){
      tempnk=sum(sum(C==k));
      tempSigma=inv(tempnk*inv(Sigma.slice(k))+inv(lambda0));
      tempMu=tempSigma*(inv(Sigma.slice(k))*sum(P(find(C==k)))+inv(lambda0)*mu0);
      Mu.col(k)=rmvnorm(1,tempMu,tempSigma);
    }
    
    //update Sigma
    for(k=0;k<K;k++){
      tempnk=sum(sum(C==k));
      tempv0=v0+tempnk;
      tempSk=computeSk(P,C,Mu,k);
      tempSk=inv(sigma0)+tempSk;
      Sigma.slice(k)=riwish(tempv0,tempSk);
    }
    
    //update Cij
    vec probK(K,fill::zeros);
    vec fullvec = regspace<vec>(1,K);
    for(i=0;i<m;i++)
      for(j=0;j<n;j++){
        for(k=0;k<K;k++){
          probK.at(k)=sum(dmvnorm(P.tube(i,j),Mu.col(k),Sigma.slice(k)))*Pr(C,alpha,beta,i,j,k);
        }
        probK/=sum(probK);
        C.at(i,j)=sum(Rcpp::RcppArmadillo::sample(fullvec,1,true,probK));
      }
    //record
    if(l>=burnIn) sampleC.slice(l-burnIn)=C;
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
