set.seed(123)
library(mvtnorm)

prediction_monte_carlo_loop = function(Mu, S, n_mc = 1e6, n_vec = seq(20, 200, 10)) {

  p = length(Mu)
  res = array(NA_real_, 
              dim=c(length(n_vec), n_mc, 7), 
              dimnames=list(
                  paste(n_vec), 
                  NULL, 
                  c('obs','m_eq','s_eq','dof_eq','m_uneq', 's_uneq', 'dof_uneq')
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
    
      # fit unequal weights
      X = cbind(1, data_tr[, -1])
      XtXm1_uneq = solve(crossprod(X))
      beta_uneq = XtXm1_uneq %*% crossprod(X, y)
      rss_uneq = sum((y - X %*% beta_uneq)^2)
      r_uneq = ncol(X) 
    
      # fit equal weights
      X = cbind(1, rowMeans(data_tr[, -1]))
      XtXm1_eq = solve(crossprod(X))
      beta_eq = XtXm1_eq %*% crossprod(X, y)
      rss_eq = sum((y - X %*% beta_eq)^2)
      r_eq = ncol(X) 
    
      # simulate test data
      data_te = drop(rmvnorm(1, mean=Mu, sigma=S))
    
      # predict unequal weights
      X_te = c(1, data_te[-1])
      m_uneq = sum(beta_uneq * X_te)
      s2_uneq = rss_uneq / (n_ - r_uneq) * (sum(tcrossprod(X_te) * XtXm1_uneq) + 1)
    
      X_te = c(1, mean(data_te[-1]))
      m_eq = sum(beta_eq * X_te)
      s2_eq = rss_eq / (n_ - r_eq) * (sum(tcrossprod(X_te) * XtXm1_eq) + 1)
       
    
      res[i_n, i_mc, ] = 
        c(obs=data_te[1],
          m_eq=m_eq,
          s_eq=sqrt(s2_eq),
          dof_eq=n_-2,
          m_uneq=m_uneq,
          s_uneq=sqrt(s2_uneq),
          dof_uneq=n_-length(Mu)-1)
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
logs_t = function(y, m, s, n) {
  return(-dt((y-m)/s, n, log=TRUE) + log(s))
}


# calculate individual scores
se_enso_eq = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_eq'])^2)
se_enso_uneq = apply(res_enso, 1, function(X) (X[,'obs'] - X[,'m_uneq'])^2)
crps_enso_eq = apply(res_enso, 1, function(X) crps_t(X[,'obs'], X[,'m_eq'], X[,'s_eq'], X[,'dof_eq']))
crps_enso_uneq = apply(res_enso, 1, function(X) crps_t(X[,'obs'], X[,'m_uneq'], X[,'s_uneq'], X[,'dof_uneq']))
logs_enso_eq = apply(res_enso, 1, function(X) logs_t(X[,'obs'], X[,'m_eq'], X[,'s_eq'], X[,'dof_eq']))
logs_enso_uneq = apply(res_enso, 1, function(X) logs_t(X[,'obs'], X[,'m_uneq'], X[,'s_uneq'], X[,'dof_uneq']))

se_nao_eq = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_eq'])^2)
se_nao_uneq = apply(res_nao, 1, function(X) (X[,'obs'] - X[,'m_uneq'])^2)
crps_nao_eq = apply(res_nao, 1, function(X) crps_t(X[,'obs'], X[,'m_eq'], X[,'s_eq'], X[,'dof_eq']))
crps_nao_uneq = apply(res_nao, 1, function(X) crps_t(X[,'obs'], X[,'m_uneq'], X[,'s_uneq'], X[,'dof_uneq']))
logs_nao_eq = apply(res_nao, 1, function(X) logs_t(X[,'obs'], X[,'m_eq'], X[,'s_eq'], X[,'dof_eq']))
logs_nao_uneq = apply(res_nao, 1, function(X) logs_t(X[,'obs'], X[,'m_uneq'], X[,'s_uneq'], X[,'dof_uneq']))

save(file='n-dependence-scores.Rdata', 
  list=c('se_enso_eq', 'se_enso_uneq',
         'logs_enso_eq', 'logs_enso_uneq',
         'crps_enso_eq', 'crps_enso_uneq',
         'se_nao_eq', 'se_nao_uneq',
         'logs_nao_eq', 'logs_nao_uneq',
         'crps_nao_eq', 'crps_nao_uneq'))

# save summary statistics
n_mc = dim(logs_nao_eq)[1]
nvec = as.numeric(dimnames(logs_nao_eq)[[2]])
logs_nao = cbind(mean_eq=colMeans(logs_nao_eq), 
                 sd_eq=apply(logs_nao_eq, 2, sd),
                 mean_uneq=colMeans(logs_nao_uneq), 
                 sd_uneq=apply(logs_nao_uneq, 2, sd))
logs_enso = cbind(mean_eq=colMeans(logs_enso_eq), 
                 sd_eq=apply(logs_enso_eq, 2, sd),
                 mean_uneq=colMeans(logs_enso_uneq), 
                 sd_uneq=apply(logs_enso_uneq, 2, sd))
se_nao = cbind(mean_eq=colMeans(se_nao_eq), 
                 sd_eq=apply(se_nao_eq, 2, sd),
                 mean_uneq=colMeans(se_nao_uneq), 
                 sd_uneq=apply(se_nao_uneq, 2, sd))
se_enso = cbind(mean_eq=colMeans(se_enso_eq), 
                 sd_eq=apply(se_enso_eq, 2, sd),
                 mean_uneq=colMeans(se_enso_uneq), 
                 sd_uneq=apply(se_enso_uneq, 2, sd))
crps_nao = cbind(mean_eq=colMeans(crps_nao_eq), 
                 sd_eq=apply(crps_nao_eq, 2, sd),
                 mean_uneq=colMeans(crps_nao_uneq), 
                 sd_uneq=apply(crps_nao_uneq, 2, sd))
crps_enso = cbind(mean_eq=colMeans(crps_enso_eq), 
                 sd_eq=apply(crps_enso_eq, 2, sd),
                 mean_uneq=colMeans(crps_enso_uneq), 
                 sd_uneq=apply(crps_enso_uneq, 2, sd))

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

       

