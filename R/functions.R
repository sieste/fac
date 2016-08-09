#' Calculate factor analysis log-likelihood function
#'

fa_llik = function(x, Mu=rep(0,p), Sigma=rep(1,p), Lambda=matrix(0, p, 1)) {
  p = ncol(x)
  Psi = 1 - rowSums(Lambda^2)
  C = (tcrossprod(Lambda) + diag(Psi)) * tcrossprod(sqrt(Sigma))
  dec = chol(C)
  tmp = backsolve(dec, t(x) - Mu, transpose=TRUE)
  rss = colSums(tmp^2)
  ans = -sum(log(diag(dec))) - .5 * p * log(2. * pi) - .5 * rss
  return(ans)
}


#' mle factor analysis a la mardia/kent/bibby sec 9.4 
#'
#' (adopted from stats:::factanal.mle.fit)
#'
fa_mle = function(x, nf=1) {
  p = ncol(x)
  n = nrow(x)

  # mles of mu and sigma
  Mu = colMeans(x)
  Sigma = colMeans(x^2) - Mu^2
  
  # correlation matrix
  S = cor(x)

  FAout = function(Psi, S, q) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L:q, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L:q] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    return(load)
  }
  FAfn = function(Psi, S, q) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE, only.values = TRUE)
    e     = E$values[-(1L:q)]
    return(-sum(log(e) - e) - q + nrow(S))
  }
  FAgr = function(Psi, S, q) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L:q, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L:q] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    g     = rowSums(load^2) + Psi - diag(S)
    return(g/Psi^2)
  }
  start = (1 - 0.5 * nf/p)/diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = 0.1, 
              upper = 1, control = list(fnscale = 1, parscale = rep(0.01, 
              length(start))), q = nf, S = S)
  Lambda = FAout(res$par, S, nf)
  dof    = p * (nf + 3) - 0.5 * nf * (nf - 1)
  Psi    = 1 - rowSums(Lambda^2)
  BIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda)) + log(n) * dof
  ans    = list(Mu = Mu, Sigma=Sigma, Lambda=Lambda, Psi=Psi,
                BIC=BIC, converged = res$convergence == 0, dof=dof)
  return(ans)
}
  
