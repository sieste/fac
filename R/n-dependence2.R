# sample sizes
n_vec = c(seq(10,100,10), seq(125,200,25))


###################################################
# load normal distribution parameters
load('../data/enso-nao.Rdata')
nao[, -1] = nao[, -1] / 1000
mu_nao = colMeans(nao)
sigma_nao = var(nao) 
mu_enso = colMeans(enso)
sigma_enso = var(enso) 

###################################################
# climatological prediction function
f_clim = function(x_tr, y_tr, x_te) {
  
  # prediction
  m_pred = mean(y_tr)
#  v_pred = var(y_tr) * (1 + 1 / length(y_tr))
#  dof_pred = length(y_tr) - 1
  v_pred = var(y_tr)
  dof_pred = Inf

  c(m_pred, v_pred, dof_pred)

}


###################################################
# multiple linear regression function
f_mlr = function(x_tr, y_tr, x_te) {

  # parameter estimation
  X_tr = cbind(1, x_tr)
  XtXm1 = solve(crossprod(X_tr))
  beta = XtXm1 %*% crossprod(X_tr, y_tr)
  y_hat = X_tr %*% beta
  rss = sum((y_tr - y_hat)^2)
  
  # prediction
  X_te = c(1, x_te)
  m_pred = sum(beta * X_te)
#  dof_pred = n - length(beta)
#  v_pred = rss / dof_pred * (1 + X_te %*% XtXm1 %*% X_te)
  v_pred = rss / (n - length(beta))
  dof_pred = Inf
  
  return(c(m_pred, v_pred, dof_pred))
  
}


###################################################
# simple linear regression function on the mmm
f_mmm = function(x_tr, y_tr, x_te) {
  
  # summary statistics
  X_tr = rowMeans(x_tr)
  X_te = mean(x_te)
  
  return(f_mlr(X_tr, y_tr, X_te))
}


###################################################
# maximum likelihood estimation for the one-factor model
fa_mle = function(x) {

  Mu = colMeans(x)
  Sigma = colMeans(x^2) - Mu^2
  
  # correlation matrix for scale-free factor analysis
  S = cor(x)

  FAout = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    return(load)
  }
  FAfn = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE, only.values = TRUE)
    e     = E$values[-1L]
    return(-sum(log(e) - e) - 1 + nrow(S))
  }
  FAgr = function(Psi, S) {
    Sstar = S * tcrossprod(1/sqrt(Psi))
    E     = eigen(Sstar, symmetric = TRUE)
    L     = E$vectors[, 1L, drop = FALSE]
    load  = L * rep(sqrt(pmax(E$values[1L] - 1, 0)), each=ncol(S)) * sqrt(Psi)
    g     = rowSums(load^2) + Psi - diag(S)
    return(g/Psi^2)
  }

  start = (1 - 0.5 / ncol(x)) / diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = 1e-6,
                upper = 1, S = S)
  if (res$convergence != 0) {
    warning('fa_mle has not converged')
  }
  Psi    = res$par
  Lambda = drop(FAout(Psi, S))
  if (Lambda[1L] < 0) {
    Lambda = -1 * Lambda
  }
  return(list(Lambda=sqrt(Sigma) * Lambda, Psi=Sigma*Psi))
}


###################################################
# predictive factor analysis regression function
f_fa = function(x_tr, y_tr, x_te) {

  # factor analysis
  fa = fa_mle(cbind(y_tr, x_tr))
  Psi = fa$Psi
  Lambda = fa$Lambda
  Lambda_x = Lambda[-1L]
  Psi_x = Psi[-1L]
  
  # regression parameters
  n = length(y_tr)
  m_x = colSums(x_tr) / nrow(x_tr)
  m_y = sum(y_tr) / length(y_tr)
  XtXinv = diag(1/Psi_x) - (tcrossprod(Lambda_x/Psi_x)) / 
                           (1 + sum(Lambda_x^2 / Psi_x))
  beta = Lambda[1L] / (1 + sum(Lambda_x^2 / Psi_x)) * Lambda_x / Psi_x
  y_hat = m_y - sum(beta * m_x) + x_tr %*% beta 
  rss = sum((y_tr - y_hat)^2)

  # prediction
  m_pred = m_y + sum(beta * (x_te - m_x))
  dof_pred = n - length(beta) - 1 
  v_pred = Psi[1L] + Lambda[1L]^2 / (1 + sum(Lambda_x^2/Psi_x))
  # v_pred = v_pred (1 + 1/n)
  dof_pred = Inf
  return(c(m_pred, v_pred, dof_pred))
}


###################################################
# score functions
f_sqerr = function(y_te, m_pred, v_pred, dof_pred) {
  return((y_te - m_pred)^2)
}
f_logs = function(y_te, m_pred, v_pred, dof_pred) {
  z = (y_te - m_pred) / sqrt(v_pred)
  return( -dt(x=z, df=dof_pred, log=TRUE) + 0.5 * log(v_pred) )
}
f_crps = function(y_te, m_pred, v_pred, dof_pred) {
  s_pred = sqrt(v_pred)
  z = (y_te - m_pred) / s_pred
  if (dof_pred == Inf) {
    ans = s_pred * (z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi))
  }  else {
    ans = s_pred * 
      (2 * dt(z, dof_pred) * (dof_pred + z^2) / (dof_pred-1) + 
      z * (2*pt(z, dof_pred) - 1) - 
      2*sqrt(dof_pred) / (dof_pred - 1) * beta(0.5, dof_pred - 0.5) / 
      beta(0.5, 0.5*dof_pred)^2)
  }
  return(ans)
}


###################################################
# simulation function
sim = function(mu, sigma, n, seed) {
  set.seed(seed)
  return(mvtnorm::rmvnorm(n, mu, sigma))
}


###################################################
###################################################
###################################################
###################################################
# monte carlo experiment

n_mc = 5e5

n_vec = c(seq(10,100,10), seq(125,200,25))
parameters = list(
  nao=list(mu=mu_nao, sigma=sigma_nao), 
  enso=list(mu=mu_enso, sigma=sigma_enso)
)
methods = list(
  clim = f_clim,
  mmm = f_mmm,
  mlr = f_mlr,
  fa = f_fa
)
scores = list(
  crps=f_crps,
  logs=f_logs,
  sqerr=f_sqerr
)

res = array(NA_real_, 
  dim = c(length(parameters),
          length(methods),
          length(scores),
          length(n_vec),
          n_mc),
  dimnames = list(names(parameters),
                  names(methods),
                  names(scores),
                  paste(n_vec),
                  NULL)
)
 
for (i_mc in 1:n_mc) {

  for (i_n in seq_along(n_vec)) {
    n = n_vec[i_n]
    nam_n = paste(n)
  
    for (i_par in seq_along(parameters)) {
  
      nam_parameter = names(parameters)[i_par]
      mu = parameters[[i_par]]$mu
      sigma = parameters[[i_par]]$sigma
  
      # simulate data, last row is test data
      data = sim(mu, sigma, n + 1, seed=i_mc)
  
      for (i_method in seq_along(methods)) {
        
        nam_method = names(methods)[i_method]
        f = methods[[i_method]]
        pred = f(x_tr=data[1:n, -1], y_tr=data[1:n, 1], x_te=data[n+1, -1])
  
        for (i_score in seq_along(scores)) {
          nam_score = names(scores)[[i_score]] 
          f_score = scores[[i_score]]
          score = f_score(data[n + 1, 1], pred[1], pred[2], pred[3])


          if (i_mc %% 1000 == 0) {
            print(c(nam_parameter, nam_method, nam_score, nam_n, i_mc))
          }
          res[nam_parameter, nam_method, nam_score, nam_n, i_mc] = score
        }
      }
    }
  }
}

save(file='n-dependence2.Rdata', list='res')

