#Sys.setenv("PKG_LIBS"="-ltrng4")
Sys.setenv("PKG_LIBS"="-ltrng4 -lprofiler")
Rcpp::sourceCpp('~/Rstudio_ccy/Imagesegment/ImageGibbs.cpp')
library(jpeg)
P = readJPEG("./seg.jpg")
#P=P[1:150,1:150,]
P_seg = P

K = 14 #number of clusters
R = P[, , 1]
G = P[, , 2]
B = P[, , 3]
m = nrow(R)
n = ncol(R)

dim(R) = c(m * n, 1)
dim(G) = c(m * n, 1)
dim(B) = c(m * n, 1)

data = cbind(R, G, B)

Kmeans = kmeans(data, K) #kmeans K clusters

# Initial value of C
C = Kmeans$cluster
dim(C) = c(m, n)

for (i in 1:m) {
  for (j in 1:n) {
    k = C[i, j]
    P_seg[i, j, ] = Kmeans$centers[k, ]
  }
}

# Initial value of Sigma,Mu
Sigma = array(1:(9 * K), dim = c(3, 3, K))
Mu = matrix(1:3 * K, 3, K)
for (k in 1:K) {
  C1 = data[Kmeans$cluster == k, ]
  Sigma[, , k] = cov(C1)
  Mu[, k] = mean(C1)
}

# Initial values
alpha = 0.8
beta = 0.2
v0 = 7 # df of prior Inverse-Wishart
sigma0 = diag(rep(1, 3))
mu0 = c(1, 1, 1)
lambda0 = diag(rep(1, 3))
Prob <- matrix(rep(0,m*n),nrow = m,ncol = n)

# library(microbenchmark)
# microbenchmark(ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 0,1),times = 1)
Result=ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 500,500)
save.image(paste("~/Rstudio_ccy/script",date(),".RData"))