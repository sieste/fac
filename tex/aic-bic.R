sqrt0 = function(x) {
  x[x<0] = 0
  sqrt(x)
}

ll = function(pars, X, Rmat, M) {

  mu_x = pars[1:M]
  sigma_x = pars[(M+1):(2*M)]
  Lambda_x = matrix(pars[(2*M+1):length(pars)], nrow=M)
  
  delta_x = sqrt0(1 - rowSums(Lambda_x^2))
  Lambda = Rmat %*% Lambda_x
  delta = drop(Rmat %*% delta_x)
  sigma = drop(Rmat %*% sigma_x)
  mu = drop(Rmat %*% mu_x)

  # multivariate normal
  Mu = mu
  Sigma = diag(sigma) %*% (tcrossprod(Lambda) + diag(delta)^2) %*% diag(sigma)

  return(-sum(mvtnorm::dmvnorm(X, Mu, Sigma, log=TRUE)))
}

load("~/folders/nino34-combination/nino34.Rdata")
M = length(lst)
Rvec = R
Rmat = matrix(0, nrow=sum(Rvec), ncol=length(Rvec))
k = 0
for (r in 1:length(Rvec)) {
  for (rr in 1:Rvec[r]) {
    k = k+1
    Rmat[k, r] = 1
  }
}


opts = list()

for (S in 1:10) {
  pars = c(
    sapply(lst, mean),
    sapply(lst, function(x) mean(apply(x, 2, sd))),
    rep(.7/sqrt(S), M*S)
  )
  ans = optim(pars, ll, method='BFGS', X=t(ens), 
              Rmat=Rmat, M=M, 
              control=list(maxit=1e5, trace=10))
  opts[[paste(S)]] = ans
  print(ans)
                     
}


k = sapply(sapply(opts, `[`, 'par'), length)
n = nrow(ens)
bic = 2* unlist(sapply(opts, `[`, 'value')) + k*log(n)
aic = 2* unlist(sapply(opts, `[`, 'value')) + k*2

save(file='aic-bic.Rdata', list=c('aic', 'bic'))
