library(jpeg)
P = readJPEG("./seg.jpg")
P_seg = P

K = 14 #number of clusters
R=P[,,1]
G=P[,,2]
B=P[,,3]
m = nrow(R)
n = ncol(R)

dim(R) = c(m*n,1)
dim(G) = c(m*n,1)
dim(B) = c(m*n,1)

data = cbind(R,G,B)

Kmeans = kmeans(data,K) #kmeans K clusters

# Initial value of C
C = Kmeans$cluster
dim(C) = c(m,n) 

for(i in 1:m){
  for(j in 1:n){
    k = C[i,j]
    P_seg[i,j,] = Kmeans$centers[k,]
  }
}

# Initial value of Sigmak,Muk
SigmaK= array(1:(9*K),dim = c(3,3,K))
MuK=matrix(1:3*K,3,K)
for(k in 1:K){
  C1 = data[Kmeans$cluster==k,]
  SigmaK[,,k]=cov(C1)
  MuK[,k]=mean(C1)
}

# Initial values
alpha=2
beta=5
v0=7 # df of prior Inverse-Wishart
sigma0=matrix(1:9,3,3)
mu0=c(1,1,1)
lambda0=matrix(rep(1,9),3,3)

result=ImageGibbs(K,P,C,Muk,SigmaK,alpha,beta,mu0,lambda0,v0,sigma0,1000,10000)

writeJPEG(P_seg,"./my1.jpg",0.95)
