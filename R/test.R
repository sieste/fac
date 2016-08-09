source('functions.R')
library(mvtnorm)
set.seed(12345)

# artificial data
n = 100
p = 10
Mu = rep(2, p)
Lambda = matrix(rnorm(2*p)^2, p, 2)
Psi = rep(1, p)
Sigma = tcrossprod(Lambda) + diag(c(Psi^2))

x = rmvnorm(n=n, mean=Mu, sigma=Sigma)

# test factor analysis log likelihood function
ll1 = fa_llik(x, Mu, Lambda, Psi)
ll2 = mvtnorm::dmvnorm(x, Mu, Sigma, log=TRUE)
stopifnot(all(ll1 == ll2))

# microbenchmark::microbenchmark(
# fa_llik(x, Mu, Lambda, Psi),
# mvtnorm::dmvnorm(x, Mu, tcrossprod(Lambda)+diag(c(Psi)), log=TRUE)
# )
# 75 vs 170 microseconds

print(fa_mle(x, control=list(maxit=1e6)))

for (nf in 1:4) {
  opt = fa_mle(x, nf)
  print(fa_bic(opt))
}

