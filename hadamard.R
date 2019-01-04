rm(list=ls())
library(microbenchmark)
library(measError)
# Parameters
n = 2^17 # number of samples
d = 50 # number of variables
m = 500
q = 1e-4
q*n

# Generate the design matrix
A = matrix(rnorm(n*d), ncol = d)
# Generate the regression coefficients
beta = rnorm(d, sd = 5)
## Generate the response
b = A %*% beta + rnorm(n, sd = 1)
Ab = cbind(A, b)

## Direct least square regression
betaHat = solve(t(A) %*% A, t(A) %*% b)
betaHat2 = randomizedRegression(Ab, m, q)
cbind(betaHat, betaHat2)

## Profiling
Sys.setenv(CPUPROFILE_FREQUENCY = 1000)
start_profiler("/home/chencanyi2018/Rstudio_ccy/profile.out")
for(i in 1:5){
  betaHat2 = randomizedRegression(Ab, m, q)
}
stop_profiler()