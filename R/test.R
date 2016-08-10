source('functions.R')
library(mvtnorm)
set.seed(12345)

# artificial data
n = 30
p = 6
nf = 3
Mu = rep(0, p)
Sigma = rep(1, p)
Lambda = gtools::rdirichlet(p, rep(1,nf))
Psi = 1 - rowSums(Lambda^2)
C = (tcrossprod(Lambda) + diag(Psi)) * tcrossprod(Sigma)

x = rmvnorm(n=n, mean=Mu, sigma=C)

# test factor analysis log likelihood function
ll1 = fa_llik(x, Mu, Sigma, Lambda)
ll2 = mvtnorm::dmvnorm(x, Mu, C, log=TRUE)
stopifnot(all(ll1 == ll2))

# microbenchmark::microbenchmark(
# fa_llik(x, Mu, Lambda, Psi),
# mvtnorm::dmvnorm(x, Mu, tcrossprod(Lambda)+diag(c(Psi)), log=TRUE)
# )
# 75 vs 170 microseconds

print(fa_mle(x, 1))
print(fa_mle(x, 2))
print(fa_mle(x, 3))

