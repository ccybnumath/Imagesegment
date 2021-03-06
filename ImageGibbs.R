library(jpeg)
P = readJPEG("./seg.jpg")
#P = readJPEG("./seg_big.jpg")
#P=P[1:150,1:150,]
P_seg = P

K = 12 #number of clusters
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


#result = ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 500, 6000)
result = ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 10, 100)



Mod <- function(x) {
  as.numeric(names(table(x)))[which.max(table(x))]
}

findMod <- function(result) {
  C <- result[, , 1]
  for (i in 1:nrow(C)) {
    for (j in 1:ncol(C)) {
      C[i, j] = Mod(result[i, j, ])
    }
  }
  return(C)
}

library(RColorBrewer)
C1 <- findMod(Result)
library(RColorBrewer)
color <- col2rgb(brewer.pal(12, "Paired"))

C1 <- C1+1
for (i in 1:m) {
  for (j in 1:n) {
    k = C1[i, j]
    k = min(k, 12)
    P_seg[i, j, ] = as.vector(color[, k])
  }
}
writeJPEG(P_seg, "./my1.jpg", 0.95)
writeJPEG(P, "./my2.jpg", 0.95)

# test speed
library(microbenchmark)
microbenchmark(computeSk_parallel(P, C, Mu, 10), computeSk(P, C, Mu, 10))
microbenchmark(sumP(P, C, 10), sumP_parallel(P, C, 10))
microbenchmark(
  UpdateC(C, P, Mu, Sigma, m, n, K, alpha, beta, 15, 19),
  UpdateC_parallel(C, P, Mu, Sigma, m, n, K, alpha, beta, 15, 19)
)
library(microbenchmark)
microbenchmark(UpdateCij(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19),UpdateCij_parallel(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19),UpdateCij_parallelUnique(C,P,Mu,Sigma,m,n,K,alpha,beta,15,19,Prob))
microbenchmark(UpdateC_parallel(C,P,Mu,Sigma,m,n,K,alpha,beta,Prob),times = 1)
microbenchmark(ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 0,1),times = 1)

#profile
Sys.setenv(CPUPROFILE_FREQUENCY = 1000)
start_profiler("/home/chencanyi2018/Rstudio_ccy/profile.out")
ImageGibbs(K, P, C, Mu, Sigma, alpha, beta, mu0, lambda0, v0, sigma0, 0,1)
stop_profiler()

pprof --text /home/chencanyi2018/Rstudio_ccy/profile.out
pprof --svg /home/chencanyi2018/Rstudio_ccy/profile.out >