source('functions.R')
set.seed(123)

# artificial data
x = matrix(rnorm(100), 20, 5)
n = nrow(x)
p = ncol(x)
Lambda = rep(.5, p)
Psi = rep(1, p)
Mu = rep(0, p)

# test factor analysis log likelihood function
ll1 = fa_llik(x, Mu, Lambda, Psi)
ll2 = mvtnorm::dmvnorm(x, Mu, tcrossprod(Lambda)+diag(c(Psi)), log=TRUE)
stopifnot(all(ll1 == ll2))

# microbenchmark::microbenchmark(
# fa_llik(x, Mu, Lambda, Psi),
# mvtnorm::dmvnorm(x, Mu, tcrossprod(Lambda)+diag(c(Psi)), log=TRUE)
# )
# 75 vs 170 microseconds



