set.seed(123)
library(mvtnorm)

fa_mle = function(x, nf=1L, Psi.lower=1e-6) {

  p = ncol(x)

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

  start = (1 - 0.5 / p) / diag(solve(S))
  res   = optim(start, FAfn, FAgr, method = "L-BFGS-B", lower = Psi.lower,
                upper = 1, S = S)

  Lambda = drop(FAout(res$par, S))
  if (Lambda[1L] < 0) {
    Lambda = -1 * Lambda
  }
  Psi    = res$par
  ans    = list(Lambda=sqrt(Sigma) * Lambda, Psi=Sigma*Psi,
                converged = res$convergence == 0)

  return(ans)
}


prediction_monte_carlo_loop = 
function(Mu, S, 
         n_mc = 1e6, 
         n_vec = c(seq(10,100,10), seq(125,200,25))) 
{

  p = length(Mu)
  res = array(NA_real_, 
              dim=c(length(n_vec), n_mc, 9), 
              dimnames=list(
                  paste(n_vec), 
                  NULL, 
                  c('obs',
                    'm_eq','s_eq',
                    'm_uneq', 's_uneq', 
                    'm_cl', 's_cl', 
                    'm_fa', 's_fa')
              ))


  for (i_n in seq_along(n_vec)) {
    for (i_mc in seq_len(n_mc)) {
      n_ = n_vec[i_n]

      if (i_mc %% 1e4 == 0) {
        print(match.call())
        cat('n = ', n_, ' i_mc = ', i_mc, '\n')
      }
  
      # simulate independent sample of size n
      data_tr = rmvnorm(n_, mean=Mu, sigma=S)
      y = data_tr[,1]
      m_y = mean(y)
      m_x = colMeans(data_tr[, -1])
    
      # fit unequal weights
      X_uneq = cbind(1, data_tr[, -1])
      XtXm1_uneq = solve(crossprod(X_uneq))
      beta_uneq = XtXm1_uneq %*% crossprod(X_uneq, y)
      rss_uneq = sum((y - X_uneq %*% beta_uneq)^2)
      r_uneq = ncol(X_uneq) 
    
      # fit equal weights
      X_eq = cbind(1, rowMeans(data_tr[, -1]))
      XtXm1_eq = solve(crossprod(X_eq))
      beta_eq = XtXm1_eq %*% crossprod(X_eq, y)
      rss_eq = sum((y - X_eq %*% beta_eq)^2)
      r_eq = ncol(X_eq)

      # factor analysis model
      fa = fa_mle(data_tr)
      beta_fa = fa$Lambda[1] / (1 + sum(fa$Lambda^2 / fa$Psi)) * fa$Lambda[-1]/fa$Psi[-1]
      y_hat = m_y - sum(m_x * beta_fa) + data_tr[, -1] %*% beta_fa
      rss_fa = sum((y - y_hat)^2)
      r_fa = ncol(data_tr)
    
      # simulate test data
      data_te = drop(rmvnorm(1, mean=Mu, sigma=S))
    
      # predict with unequal weights
      X_te = c(1, data_te[-1])
      m_uneq = sum(beta_uneq * X_te)
      s2_uneq = rss_uneq / (n_ - r_uneq) 
    
      # predict with equal weights
      X_te = c(1, mean(data_te[-1]))
      m_eq = sum(beta_eq * X_te)
      s2_eq = rss_eq / (n_ - r_eq) 

      # predict with factor analysis
      X_te = data_te[-1]
      m_fa = m_y - sum(m_x * beta_fa) + sum(data_te[-1] * beta_fa)
      s2_fa = rss_fa / (n_ - r_fa) 
       
      # predict with climatology
      m_cl = mean(y)
      s2_cl = var(y)
    
      res[i_n, i_mc, ] = 
        c(obs=data_te[1],
          m_eq=m_eq,
          s_eq=sqrt(s2_eq),
          m_uneq=m_uneq,
          s_uneq=sqrt(s2_uneq),
          m_cl=m_cl,
          s_cl=sqrt(s2_cl),
          m_fa=m_fa,
          s_fa=sqrt(s2_fa))
    }
  }

  return(res)

}

# enso and nao data
load('../data/enso-nao.Rdata')

# ENSO analysis
n_enso = nrow(enso)
S_enso = var(enso) * (n_enso-1) / n_enso
mu_enso = colMeans(enso)

res_enso = prediction_monte_carlo_loop(mu_enso, S_enso)

# NAO analysis
n_nao = nrow(nao)
S_nao = var(nao) * (n_nao-1) / n_nao
mu_nao = colMeans(nao)

res_nao = prediction_monte_carlo_loop(mu_nao, S_nao)

save(file='n-dependence.Rdata', list=c('res_enso', 'res_nao'))




# CRPS of the non-standardised t-distribution (jordan 2015
# https://github.com/FK83/scoringRules/blob/master/crps.pdf)
crps_t = function(y, m, s, n) { 
  z = (y - m) / s
  return(s * (2 * dt(z, n) * (n + z^2) / (n-1) + z * (2 * pt(z, n) - 1) - 2 * sqrt(n) / (n-1) * beta(.5, n-.5) / beta(.5, .5*n)^2))
}
crps_n = function(y, m, s) { 
  z = (y - m) / s
  return(s * (2 * dnorm(z) + z * (2 * pnorm(z) - 1) - 1 / sqrt(pi)))
}
logs_t = function(y, m, s, n) {
  return(-dt((y-m)/s, n, log=TRUE) + log(s))
}
logs_n = function(y, m, s) {
  return(-dnorm((y-m)/s, log=TRUE) + log(s))
}


# calculate individual scores
se_enso_eq = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_eq'])^2)
se_enso_uneq = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_uneq'])^2)
se_enso_cl = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_cl'])^2)
se_enso_fa = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_fa'])^2)

crps_enso_eq = apply(res_enso, 1, function(X) crps_n(X[,'obs'], X[,'m_eq'], X[,'s_eq']))
crps_enso_uneq = apply(res_enso, 1, function(X) crps_n(X[,'obs'], X[,'m_uneq'], X[,'s_uneq']))
crps_enso_cl = apply(res_enso, 1, function(X) crps_n(X[,'obs'], X[,'m_cl'], X[,'s_cl']))
crps_enso_fa = apply(res_enso, 1, function(X) crps_n(X[,'obs'], X[,'m_fa'], X[,'s_fa']))

logs_enso_eq = apply(res_enso, 1, function(X) logs_n(X[,'obs'], X[,'m_eq'], X[,'s_eq']))
logs_enso_uneq = apply(res_enso, 1, function(X) logs_n(X[,'obs'], X[,'m_uneq'], X[,'s_uneq']))
logs_enso_cl = apply(res_enso, 1, function(X) logs_n(X[,'obs'], X[,'m_cl'], X[,'s_cl']))
logs_enso_fa = apply(res_enso, 1, function(X) logs_n(X[,'obs'], X[,'m_fa'], X[,'s_fa']))

se_nao_eq = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_eq'])^2)
se_nao_uneq = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_uneq'])^2)
se_nao_cl = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_cl'])^2)
se_nao_fa = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_fa'])^2)

crps_nao_eq = apply(res_nao, 1, function(X) crps_n(X[,'obs'], X[,'m_eq'], X[,'s_eq']))
crps_nao_uneq = apply(res_nao, 1, function(X) crps_n(X[,'obs'], X[,'m_uneq'], X[,'s_uneq']))
crps_nao_cl = apply(res_nao, 1, function(X) crps_n(X[,'obs'], X[,'m_cl'], X[,'s_cl']))
crps_nao_fa = apply(res_nao, 1, function(X) crps_n(X[,'obs'], X[,'m_fa'], X[,'s_fa']))

logs_nao_eq = apply(res_nao, 1, function(X) logs_n(X[,'obs'], X[,'m_eq'], X[,'s_eq']))
logs_nao_uneq = apply(res_nao, 1, function(X) logs_n(X[,'obs'], X[,'m_uneq'], X[,'s_uneq']))
logs_nao_cl = apply(res_nao, 1, function(X) logs_n(X[,'obs'], X[,'m_cl'], X[,'s_cl']))
logs_nao_fa = apply(res_nao, 1, function(X) logs_n(X[,'obs'], X[,'m_fa'], X[,'s_fa']))


save(file='n-dependence-scores.Rdata', 
  list=c('se_enso_eq', 'se_enso_uneq', 'se_enso_cl', 'se_enso_fa',
         'logs_enso_eq', 'logs_enso_uneq', 'logs_enso_cl', 'logs_enso_fa',
         'crps_enso_eq', 'crps_enso_uneq', 'crps_enso_cl', 'crps_enso_fa',
         'se_nao_eq', 'se_nao_uneq', 'se_nao_cl', 'se_nao_fa',
         'logs_nao_eq', 'logs_nao_uneq', 'logs_nao_cl', 'logs_nao_fa',
         'crps_nao_eq', 'crps_nao_uneq', 'crps_nao_cl', 'crps_nao_fa'))

# save summary statistics
n_mc = dim(logs_nao_eq)[1]
nvec = as.numeric(dimnames(logs_nao_eq)[[2]])
logs_nao = cbind(mean_eq=colMeans(logs_nao_eq), 
                 sd_eq=apply(logs_nao_eq, 2, sd),
                 mean_uneq=colMeans(logs_nao_uneq), 
                 sd_uneq=apply(logs_nao_uneq, 2, sd),
                 mean_fa=colMeans(logs_nao_fa), 
                 sd_fa=apply(logs_nao_fa, 2, sd),
                 mean_cl=colMeans(logs_nao_cl), 
                 sd_cl=apply(logs_nao_cl, 2, sd))
logs_enso = cbind(mean_eq=colMeans(logs_enso_eq), 
                 sd_eq=apply(logs_enso_eq, 2, sd),
                 mean_uneq=colMeans(logs_enso_uneq), 
                 sd_uneq=apply(logs_enso_uneq, 2, sd),
                 mean_fa=colMeans(logs_enso_fa), 
                 sd_fa=apply(logs_enso_fa, 2, sd),
                 mean_cl=colMeans(logs_enso_cl), 
                 sd_cl=apply(logs_enso_cl, 2, sd))
se_nao = cbind(mean_eq=colMeans(se_nao_eq), 
                 sd_eq=apply(se_nao_eq, 2, sd),
                 mean_uneq=colMeans(se_nao_uneq), 
                 sd_uneq=apply(se_nao_uneq, 2, sd),
                 mean_fa=colMeans(se_nao_fa), 
                 sd_fa=apply(se_nao_fa, 2, sd),
                 mean_cl=colMeans(se_nao_cl), 
                 sd_cl=apply(se_nao_cl, 2, sd))
se_enso = cbind(mean_eq=colMeans(se_enso_eq), 
                 sd_eq=apply(se_enso_eq, 2, sd),
                 mean_uneq=colMeans(se_enso_uneq), 
                 sd_uneq=apply(se_enso_uneq, 2, sd),
                 mean_fa=colMeans(se_enso_fa), 
                 sd_fa=apply(se_enso_fa, 2, sd),
                 mean_cl=colMeans(se_enso_cl), 
                 sd_cl=apply(se_enso_cl, 2, sd))
crps_nao = cbind(mean_eq=colMeans(crps_nao_eq), 
                 sd_eq=apply(crps_nao_eq, 2, sd),
                 mean_uneq=colMeans(crps_nao_uneq), 
                 sd_uneq=apply(crps_nao_uneq, 2, sd),
                 mean_fa=colMeans(crps_nao_fa), 
                 sd_fa=apply(crps_nao_fa, 2, sd),
                 mean_cl=colMeans(crps_nao_cl), 
                 sd_cl=apply(crps_nao_cl, 2, sd))
crps_enso = cbind(mean_eq=colMeans(crps_enso_eq), 
                 sd_eq=apply(crps_enso_eq, 2, sd),
                 mean_uneq=colMeans(crps_enso_uneq), 
                 sd_uneq=apply(crps_enso_uneq, 2, sd),
                 mean_fa=colMeans(crps_enso_fa), 
                 sd_fa=apply(crps_enso_fa, 2, sd),
                 mean_cl=colMeans(crps_enso_cl), 
                 sd_cl=apply(crps_enso_cl, 2, sd))

impr_nao = cbind(se=colMeans(se_nao_uneq < se_nao_eq),
                 crps=colMeans(crps_nao_uneq < crps_nao_eq),
                 logs=colMeans(logs_nao_uneq < logs_nao_eq))
impr_enso = cbind(se=colMeans(se_enso_uneq < se_enso_eq),
                 crps=colMeans(crps_enso_uneq < crps_enso_eq),
                 logs=colMeans(logs_enso_uneq < logs_enso_eq))

save(file='n-dependence-scores-summary.Rdata',
  list=c('logs_nao', 'se_nao', 'crps_nao',
         'logs_enso', 'se_enso', 'crps_enso', 
         'impr_nao', 'impr_enso',
         'n_mc', 'nvec'))

       

