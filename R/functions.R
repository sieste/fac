#' Calculate factor analysis log-likelihood function
#'

fa_llik = function(x, Mu=rep(0,p), Sigma=rep(1,p), Lambda=matrix(0, p, 1), Psi=NULL) {
  p = ncol(x)
  if (is.null(Psi)) {
    Psi = 1 - rowSums(Lambda^2)
  }
  C = (tcrossprod(Lambda) + diag(Psi)) * tcrossprod(sqrt(Sigma))
  dec = chol(C)
  tmp = backsolve(dec, t(x) - Mu, transpose=TRUE)
  rss = colSums(tmp^2)
  ans = -sum(log(diag(dec))) - .5 * p * log(2. * pi) - .5 * rss
  return(ans)
}


#' mle factor analysis a la mardia/kent/bibby sec 9.4 
#'
#' (adopted from stats:::factanal.mle.fit, 30% faster by optimising matrix
#'  operations)
#'
fa_mle = function(x, nf=1, scaled=FALSE, Psi.fixed=NULL) {

  p = ncol(x)
  n = nrow(x)

  # mles of mu and sigma
  if (scaled) {
    x = scale(x)
    Mu = rep(0, p)
    Sigma = rep(1, p)
  } else {
    Mu = colMeans(x)
    Sigma = colMeans(x^2) - Mu^2
  }
  
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

  if (missing(Psi.fixed)) {
    start = (1 - 0.5 * nf/p)/diag(solve(S))
    res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = 0.1, 
                  upper = 1, control = list(fnscale = 1, parscale = rep(0.01, 
                  length(start))), q = nf, S = S)
  } else {
    res = list(par=Psi.fixed, convergence=0)
  }

  Lambda = FAout(res$par, S, nf)
  dof    = length(Lambda) - 0.5 * nf * (nf - 1) +
           ifelse(scaled, 0, length(Mu)+length(Sigma)) + 
           ifelse(missing(Psi.fixed), 0, -length(Psi.fixed))
  Psi    = res$par
  AIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + 2. * dof
  BIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + log(n) * dof
  ans    = list(Mu = Mu, Sigma=Sigma, Lambda=Lambda, Psi=Psi,
                AIC=AIC, BIC=BIC, converged = res$convergence == 0, dof=dof)

  return(ans)
}
 
#' reconstruct convariance matrix from factor analysis
fa_to_sigma = function(fa) {
  C = with(fa, (tcrossprod(Lambda)+diag(Psi)) * tcrossprod(sqrt(Sigma)))
  return(C)
}


#' Model selection using leave-one-out predictive density
#'
fa_loo = function(x, y, nf=1) {

  data = cbind(y, x)
  n = length(y)
  score = 0
  for (i in 1:n) {

    # covariance matrix from MLE FA
    fa = fa_mle(data[-i, ], nf=nf)
    M = fa[['Mu']]
    C = fa_to_sigma(fa)

    # regression
    M_x = M[-1]
    M_y = M[1]
    C_yy = C[1, 1, drop=TRUE]
    C_yx = C[1, -1, drop=FALSE]
    C_xx_inv = solve(C[-1, -1])
    C_xy = C[-1, 1, drop=FALSE]
    
    M_p = M_y + drop(C_yx %*% C_xx_inv %*% (x[i, ] - M_x))
    C_p = C_yy - drop(C_yx %*% C_xx_inv %*% C_xy)

    # increment leave-one-out score
    score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)
  }
  return(score/n)
}

