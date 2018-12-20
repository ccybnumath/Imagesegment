library(jpeg)
P = readJPEG("./seg.jpg")
P_seg = P

K = 5 #number of clusters
R=P[,,1]
G=P[,,2]
B=P[,,3]
m = nrow(R)
n = ncol(R)

dim(R) = c(m*n,1)
dim(G) = c(m*n,1)
dim(B) = c(m*n,1)

data = cbind(R,G,B)

Kmeans = kmeans(data,K) #kmeans 5 clusters
C = Kmeans$cluster
dim(C) = c(m,n) 

for(i in 1:m){
  for(j in 1:n){
    k = C[i,j]
    P_seg[i,j,] = Kmeans$centers[k,]
  }
}

Sigma_bar = array(1:(9*K),dim = c(3,3,K))
for(k in 1:K){
  C1 = data[Kmeans$cluster==k,]
  Sigma_bar[,,k]=cov(C1)
}


writeJPEG(P_seg,"./my1.jpg",0.95)
