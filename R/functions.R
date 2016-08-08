#' Calculate factor analysis log-likelihood function
#'

fa_llik = function(x, Mu=rep(0,p), Lambda=matrix(0, p, 1), Psi=rep(1,p)) {
  p = ncol(x)
  Sigma = tcrossprod(Lambda) + diag(c(Psi^2))
  dec = chol(Sigma)
  tmp = backsolve(dec, t(x) - Mu, transpose=TRUE)
  rss = colSums(tmp^2)
  ans = -sum(log(diag(dec))) - .5 * p * log(2. * pi) - .5 * rss
  return(ans)
}

fa_llik_opt = function(par, x, n, p, nf) {
  Mu = par[1:p]
  Lambda = matrix(par[(p+1):((nf+1)*p)], ncol=nf)
  Psi = par[((nf+1)*p+1):length(par)]
  return(-sum(fa_llik(x=x, Mu=Mu, Lambda=Lambda, Psi=Psi)))
}

fa_mle = function(x, nf=1, par=c(rep(1, (nf+1)*p), rep(1, p)), control=list()) {
  n = nrow(x)
  p = ncol(x)
  opt = optim(par=par, fn=fa_llik_opt, method='BFGS', x=x, n=n, p=p, nf=nf, control=control)
  ans = list()
  ans$Mu = opt$par[1:p]
  ans$Lambda = matrix(opt$par[(p+1):((nf+1)*p)], ncol=nf)
  ans$Psi2 = opt$par[((nf+1)*p+1):length(opt$par)]^2.
  ans$convergence = opt$convergence
  ans$nll = opt$value
  ans$n = nrow(x)
  return(ans)
}

fa_bic = function(opt) {
  nf = ncol(opt$Lambda)
  p = nrow(opt$Lambda)
  k = (nf+2) * p - nf * (nf - 1) / 2
  return(2. * opt$nll + log(opt$n) * k)
}

