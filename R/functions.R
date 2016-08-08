#' Calculate factor analysis log-likelihood function
#'

fa_llik = function(x, Mu=rep(0,p), Lambda=rep(0,p), Psi=diag(p)) {
  p = ncol(x)
  Sigma = tcrossprod(Lambda) + diag(c(Psi))
  dec = chol(Sigma)
  tmp = backsolve(dec, t(x) - Mu, transpose=TRUE)
  rss = colSums(tmp^2)
  ans = -sum(log(diag(dec))) - .5 * p * log(2. * pi) - .5 * rss
  return(ans)
}


