# init
set.seed(123)

# define number of ensemble members in the MME and corresponding R matrix
Rvec = c(1,2,3)
Rmat = matrix(0, nrow=sum(Rvec), ncol=length(Rvec))
k = 0
for (r in 1:length(Rvec)) {
  for (rr in 1:Rvec[r]) {
    k = k+1
    Rmat[k, r] = 1
  }
}

# framework specifications (model-level)
S = 1
M = length(Rvec)
Lambda_x = matrix(NA_real_, nrow=M, ncol=S)
for (m in 1:M) {
  Lambda_x[m, ] = sqrt(gtools::rdirichlet(S+1, rep(1, S+1))[1:S])
}
delta_x = sqrt(1 - rowSums(Lambda_x^2))
sigma_x = runif(M, 1, 3)
mu_x = rnorm(M)

# framework specifications (member-level)
Lambda = Rmat %*% Lambda_x
delta = drop(Rmat %*% delta_x)
sigma = drop(Rmat %*% sigma_x)
mu = drop(Rmat %*% mu_x)

# multivariate normal
Mu = mu
Sigma = diag(sigma) %*% (tcrossprod(Lambda) + diag(delta)^2) %*% diag(sigma)

# simulate data
N = 10
X = mvtnorm::rmvnorm(N, Mu, Sigma)

# check the inverse of Sigma
Sigma_inv1 = solve(Sigma)

LLt = tcrossprod(Lambda_x)
D = tcrossprod(1/sigma_x/delta_x^2)
Mmat = - D * solve(diag(M) + LLt %*% diag(Rvec/delta_x^2)) %*% LLt
tau_x = 1/delta_x^2 / sigma_x^2
Sigma_inv2 = Rmat %*% Mmat %*% t(Rmat) + diag(c(Rmat %*% tau_x))


stopifnot(max(Sigma_inv1 - Sigma_inv2) < 1e-10)
# difference is of order 1e-15

# check timing, direct inversion using solve() is faster for few small ensembles, but slower for many and/or large ensemble sizes
microbenchmark::microbenchmark(
Sigma_inv1 = solve(Sigma),
{
LLt = tcrossprod(Lambda_x)
D = tcrossprod(1/sigma_x/delta_x^2)
Mmat = - solve(diag(rep(1,M)) + LLt %*% diag(Rvec/delta_x^2)) %*% LLt * D
tau_x = 1/delta_x^2 / sigma_x^2;
Sigma_inv2 = Rmat %*% Mmat %*% t(Rmat) + diag(rep(tau_x, Rvec));}
)

# calculate tr(Sigma^{-1} S)
Smat = crossprod(X) / N - tcrossprod(colMeans(X))

trSigmaInvS_1 = sum(solve(Sigma) * Smat)

RSR = t(Rmat) %*% Smat %*% Rmat 
RdiagSR = tapply(diag(Smat), rep(1:M, Rvec), sum) 
trSigmaInvS_2 = sum(Mmat * RSR) + sum(RdiagSR * tau_x)

stopifnot(trSigmaInvS_1 - trSigmaInvS_2 < 1e-10)

# calculate log|Sigma|
logdetsigma_1 = log(det(Sigma))
logdetsigma_2 = as.numeric(determinant(diag(rep(1,S)) + crossprod(Lambda_x * sqrt(Rvec)/delta_x))$modulus) + sum(2*Rvec*log(sigma_x * delta_x)) 
stopifnot(abs(logdetsigma_1 - logdetsigma_2) < 1e-14)
microbenchmark::microbenchmark(
  determinant(Sigma)$modulus,
  determinant(diag(rep(1,S)) + crossprod(Lambda_x * sqrt(Rvec) / delta_x))$modulus + sum(2*Rvec*log(sigma_x * delta_x)) 
)

# again, the timing of the direct calculation using log(det()) is better for
# small ensemble sizes, but the transformation gets faster for larger ensemble
# sizes.

# (x-mu)'*sigma^-1*(x-mu)
xm = colMeans(X)
xm_x = tapply(xm, rep(1:M, Rvec), sum)
xm2_x = tapply(xm^2, rep(1:M, Rvec), sum)

xmusigma_1 = drop(t(xm - mu) %*% Sigma_inv1 %*% (xm - mu))
xmusigma_2 = drop(t(xm_x - Rvec*mu_x) %*% Mmat %*% (xm_x - Rvec*mu_x)) + sum(tau_x * xm2_x - 2 * tau_x * mu_x * xm_x + Rvec * tau_x * mu_x^2)

# full likelihood calculation
ll_1 = sum(mvtnorm::dmvnorm(X, Mu, Sigma, log=TRUE))
ll_2 = -N/2 * (sum(Rvec) * log(2*pi) + logdetsigma_2 + trSigmaInvS_2 + xmusigma_2)

stopifnot(ll_1 - ll_2 < 1e-14)


# compare timings of likelihood calculations
rm(list=ls())

# define number of ensemble members in the MME and corresponding R matrix
Rvec = c(1,2,3)
Rmat = matrix(0, nrow=sum(Rvec), ncol=length(Rvec))
k = 0
for (r in 1:length(Rvec)) {
  for (rr in 1:Rvec[r]) {
    k = k+1
    Rmat[k, r] = 1
  }
}

# framework specifications (model-level)
S = 1
M = length(Rvec)
Lambda_x = matrix(NA_real_, nrow=M, ncol=S)
for (m in 1:M) {
  Lambda_x[m, ] = sqrt(gtools::rdirichlet(S+1, rep(1, S+1))[1:S])
}
delta_x = sqrt(1 - rowSums(Lambda_x^2))
sigma_x = runif(M, 1, 3)
mu_x = rnorm(M)

# framework specifications (member-level)
Lambda = Rmat %*% Lambda_x
delta = drop(Rmat %*% delta_x)
sigma = drop(Rmat %*% sigma_x)
mu = drop(Rmat %*% mu_x)

# multivariate normal
Mu = mu
Sigma = diag(sigma) %*% (tcrossprod(Lambda) + diag(delta)^2) %*% diag(sigma)

# simulate data
N = 10
X = mvtnorm::rmvnorm(N, Mu, Sigma)



ll_1 = function(X, mu_x, sigma_x, Lambda_x, Rmat) {
  
  delta_x = sqrt(1 - rowSums(Lambda_x^2))
  Lambda = Rmat %*% Lambda_x
  delta = drop(Rmat %*% delta_x)
  sigma = drop(Rmat %*% sigma_x)
  mu = drop(Rmat %*% mu_x)

  # multivariate normal
  Mu = mu
  Sigma = diag(sigma) %*% (tcrossprod(Lambda) + diag(delta)^2) %*% diag(sigma)

  sum(mvtnorm::dmvnorm(X, Mu, Sigma, log=TRUE))
}

ll_2 = function(X, mu_x, sigma_x, Lambda_x, Rmat) {
  
  delta_x = sqrt(1 - rowSums(Lambda_x^2))
  Lambda = Rmat %*% Lambda_x
  delta = drop(Rmat %*% delta_x)
  sigma = drop(Rmat %*% sigma_x)
  mu = drop(Rmat %*% mu_x)

  N = nrow(X)
  S = ncol(Lambda_x)
  M = nrow(Lambda_x)
  Rvec = colSums(Rmat)
  logdetsigma = as.numeric(determinant(diag(rep(1,S)) + crossprod(Lambda_x * sqrt(Rvec)/delta_x))$modulus) + sum(2*Rvec*log(sigma_x * delta_x)) 

  Smat = crossprod(X) / N - tcrossprod(colMeans(X))
  RSR = t(Rmat) %*% Smat %*% Rmat 
  RdiagSR = tapply(diag(Smat), rep(1:M, Rvec), sum) 
  LLt = tcrossprod(Lambda_x)
  D = tcrossprod(1/sigma_x/delta_x^2)
  Mmat = - D * solve(diag(rep(1,M)) + LLt %*% diag(Rvec/delta_x^2)) %*% LLt
  tau_x = 1/delta_x^2 / sigma_x^2
  trSigmaInvS = sum(Mmat * RSR) + sum(RdiagSR * tau_x)

  xm = colMeans(X)
  xm_x = tapply(xm, rep(1:M, Rvec), sum)
  xm2_x = tapply(xm^2, rep(1:M, Rvec), sum)
  xmusigma = drop(t(xm_x - Rvec*mu_x) %*% Mmat %*% (xm_x - Rvec*mu_x)) + sum(tau_x * xm2_x - 2 * tau_x * mu_x * xm_x + Rvec * tau_x * mu_x^2)

  return(-N/2 * (sum(Rvec) * log(2*pi) + logdetsigma + trSigmaInvS + xmusigma))
}

stopifnot(abs(ll_1(X, mu_x, sigma_x, Lambda_x, Rmat) - ll_2(X, mu_x, sigma_x, Lambda_x, Rmat)) < 1e-14)

microbenchmark::microbenchmark(
ll_1(X, mu_x, sigma_x, Lambda_x, Rmat),
ll_2(X, mu_x, sigma_x, Lambda_x, Rmat)
)


