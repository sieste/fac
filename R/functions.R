#' Calculate factor analysis log-likelihood function
#'
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
fa_mle = function(x, nf=1L, scaled=FALSE, Psi.lower=1e-6, warn.overpar=FALSE) {

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

  start = (1 - 0.5 * nf/p)/diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = Psi.lower,
                upper = 1, control = list(fnscale = 1, parscale = rep(0.01, 
                length(start))), q = nf, S = S)

  Lambda = FAout(res$par, S, nf)
  rownames(Lambda) = colnames(x)
  dof    = length(Lambda) - 0.5 * nf * (nf - 1) +
           ifelse(scaled, 0, length(Mu)+length(Sigma))
  Psi    = res$par
  AIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + 2. * dof
  BIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + log(n) * dof
  s      = .5 * (p - nf)^2 - .5 * (p + nf)
  if (warn.overpar & s < 0) {
    warning('overparametrised model')
  }
  ans    = list(Mu = Mu, Sigma=Sigma, Lambda=Lambda, Psi=Psi,
                AIC=AIC, BIC=BIC, converged = res$convergence == 0, dof=dof, s=s)

  return(ans)
}
 
#' reconstruct convariance matrix from factor analysis
fa_to_sigma = function(fa) {
  C = with(fa, (tcrossprod(Lambda)+diag(Psi)) * tcrossprod(sqrt(Sigma)))
  return(C)
}


#' Model selection for factor analysis models using leave-one-out predictive density
#'
#' @param fa_fun a function that takes as input a data matrix, and returns a 
#'               list with elements Mu (vector), Sigma (vector), Lambda
#'               (matrix), Psi (vector)
#'
fa_loo = function(x, y, fa_fun, nf=1L, ...) {

  data = cbind(y, x)
  n = length(y)
  score = 0
  for (i in 1L:n) {

    # covariance matrix from MLE FA
    fa = fa_fun(data[-i, ], nf=nf, ...)
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


#' VARIMAX rotation of factor loadings
#'
varimax = function(Lambda) {
  p = nrow(Lambda)
  nf = ncol(Lambda)
  if (nf <= 1) {
    return(Lambda)
  }

  # scale rows
  sc = sqrt(rowSums(Lambda^2))
  Lambda = Lambda / sc

  # iterative algorithm by 
  TT = diag(nf)
  d = 0
  for (i in 1L:1000L) {
    z = Lambda %*% TT
    B = crossprod(Lambda, z^3 - z * rep(rep(1/p, p) %*% z^2, each=p))
    sB = La.svd(B)
    TT = sB$u %*% sB$vt
    dpast = d
    d = sum(sB$d)
    if (d < dpast * (1+1e-5)) {
      break
    }
  }
  z = Lambda %*% TT
  z = z * sc
  dimnames(z) = dimnames(Lambda)
  return(z)
}


#' mle factor analysis with spherical noise
#' a la stoica/jansson 2008 doi:10.1016/j.sigpro.2009.01.002
#' 
fa_mle_sph = function(x, nf=1L, warn.overpar=FALSE) {

  p = ncol(x)
  n = nrow(x)

  # mles of mu and sigma
  Mu = colMeans(x)
  Sigma = colMeans(x^2) - Mu^2
  
  # correlation matrix
  R = cor(x)
  eig = eigen(R)
  l = eig$values[1L:nf]
  Psi = sum(eig$values[(nf+1L):p]) / (p - nf)
  A = eig$vectors[, 1L:nf, drop=FALSE]
  Lambda = A * rep(sqrt(l-Psi), each=p)

  Psi = rep(Psi, p)
  names(Psi) = colnames(x)
  rownames(Lambda) = colnames(x)
  dof    = length(Lambda) - 0.5 * nf * (nf - 1) + length(Mu) + length(Sigma)
  AIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + 2. * dof
  BIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + log(n) * dof
  s      = .5 * (p - nf)^2 - .5 * (p + nf)
  if (warn.overpar & s < 0) {
    warning('overparametrised model')
  }
  ans    = list(Mu=Mu, Sigma=Sigma, Lambda=Lambda, Psi=Psi,
                AIC=AIC, BIC=BIC, converged=TRUE, dof=dof, s=s)

  return(ans)
  
}


#' penalised likelihood estimation to enforce similarity of the specific variances
#'
fa_ple = function(x, nf=1L, scaled=FALSE, Psi.lower=1e-6, C=0, warn.overpar=FALSE) {
  p = ncol(x)
  n = nrow(x)

  # mles of mu and sigma
  Mu = colSums(x) / n
  Sigma = colSums(x^2) / n - Mu^2
  
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
    nll = -sum(log(e) - e) - q + p + C*var(Psi)*(p-1)/p
    return(nll)
  }
  FAgr = function(Psi, S, q) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L:q, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L:q] - 1, 0)), each=p) * sqrt(Psi)
    g     = rowSums(load^2) + Psi - diag(S)
    return(g/Psi^2 + 2*C/p*(Psi-mean(Psi)))
  }

  start = (1 - 0.5 * nf/p)/diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = Psi.lower,
                upper = 1, control = list(fnscale = 1, parscale = rep(0.01, 
                length(start))), q = nf, S = S)

  Lambda = FAout(res$par, S, nf)
  rownames(Lambda) = colnames(x)
  dof    = length(Lambda) - 0.5 * nf * (nf - 1) +
           ifelse(scaled, 0, length(Mu)+length(Sigma))
  Psi    = res$par
  AIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + 2. * dof
  BIC    = -2. * sum(fa_llik(x, Mu, Sigma, Lambda, Psi)) + log(n) * dof
  s      = .5 * (p - nf)^2 - .5 * (p + nf)
  if (warn.overpar & s < 0) {
    warning('overparametrised model')
  }
  ans    = list(Mu = Mu, Sigma=Sigma, Lambda=Lambda, Psi=Psi,
                AIC=AIC, BIC=BIC, converged = res$convergence == 0, dof=dof, s=s)

  return(ans)

}


#' 1-factor prediction function
fa_pred = function(x_tr, y_tr, x_te) {

  X_tr = cbind(y_tr, x_tr)
  p = ncol(X_tr)

  # mean and scale
  Mu = colMeans(X_tr)
  Sigma = colMeans(X_tr^2) - Mu^2
  
  # correlation matrix for scale-free factor analysis
  S = cor(X_tr)

  # function to calculate max-lik Lambda for a given Psi
  FAout = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    return(load)
  }
  # function to calculate likelihood for a given Psi
  FAfn = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE, only.values = TRUE)
    e     = E$values[-1L]
    return(-sum(log(e) - e) - 1 + nrow(S))
  }
  # gradient of likelihood function
  FAgr = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    g     = rowSums(load^2) + Psi - diag(S)
    return(g/Psi^2)
  }

  # maximise likelihood
  start = (1 - 0.5 / p) / diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = 1e-6,
                upper = 1, S = S)

  # reconstruct Lambda and Psi
  Lambda = drop(FAout(res$par, S))
  if (Lambda[1L] < 0) {
    Lambda = -1 * Lambda
  }
  Lambda = sqrt(Sigma) * Lambda
  Psi = Sigma * Psi


  # parameter estimation for post-processing
  m_x = colMeans(x_tr)
  m_y = mean(y_tr)
  beta_fa = Lambda[1] / (1 + sum(Lambda^2 / Psi)) * Lambda[-1]/Psi[-1]
  y_hat = m_y - sum(m_x * beta_fa) + x_tr %*% beta_fa
  rss_fa = sum((y - y_hat)^2)

  # prediction
  m_pred = m_y + sum(beta_fa * (x_te - m_x))
  # calculate (t(X) %*% X)^{-1} using maximum likelihood estimators of factor
  # analysis; the inflation factor is (x0' (X'X)^-1 x0 + 1)
  v_inflation = 1 + drop(x_te %*% solve(crossprod(X_tr)) %*% x_te)
  dof = length(beta_fa) + 1
  v_pred = rss_fa / (nrow(X_tr) - dof)

  # return
  return(c(m_pred, v_pred, v_inflation, dof))
}

