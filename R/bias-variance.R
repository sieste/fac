# simulate with known parameters

# estimate the parameters with different methods:

# unequal weighting, unbiased MLE estimators

# equal weighting, biased estimators

# some form of shrinkage estimator, e.g. factor analysis or principal component
# regression

# do it many times for different training sample sizes and parameter settings

# analyse bias, variance and error of parameter estimates

# check how this influences out-of-sample prediction scores


############################################################

# calculate enso and nao parameters
load('../data/enso-nao.Rdata')
n_enso = nrow(enso)
m_enso = ncol(enso) - 1
p_enso = ncol(enso)
mu_enso = colMeans(enso)
Sigma_enso = var(enso) * (n_enso - 1) / n_enso
n_nao = nrow(nao)
m_nao = ncol(nao) - 1
p_nao = ncol(nao)
mu_nao = colMeans(nao)
Sigma_nao = var(nao) * (n_nao - 1) / n_nao



sim = function(Mu, Sigma, N) {
  return(mvtnorm::rmvnorm(N, Mu, Sigma))
}

est_weights = function(data, method=c('equal', 'unequal')) {

  y = data[,1]
  if (method == 'equal') {
    X = rowMeans(data[, -1])
    weights = solve(crossprod(X)) %*% crossprod(X, y)
    weights = weights / (ncol(data) - 1)
  }
  if (method == 'unequal') {
    X = data[, -1]
    weights = solve(crossprod(X)) %*% crossprod(X, y)
  }
  return(weights)

}

# calculate the true weights mean(y|x) = mean(y) + Sigma_yx * inv(Sigma_xx) *
# (x - mean(x)); we interpret Sigma_yx * inv(Sigma_xx) as the vector of
# 'weights'
true_weights = function(Sigma) {
  Sigma[1, -1] %*% solve(Sigma[-1, -1])
}


n_vec = c(seq(10,100,10), seq(125,200,25))
n_mc = 1e5
set.seed(123)
lst_nao = list()
lst_enso = list()
for (i_n in seq_along(n_vec)) {

  n_ = n_vec[i_n]
  print(n_)

  # NAO
  w = replicate(n_mc, {
    data = sim(mu_nao, Sigma_nao, n_)
    w_eq = est_weights(data, 'equal')
    w_uneq = est_weights(data, 'unequal')
    c(w_eq, w_uneq)
  })

  w_est_eq = matrix(w[1, ], nrow=m_nao, ncol=n_mc, byrow=TRUE)
  w_est_uneq = w[-1, ]
  w_true = matrix(true_weights(Sigma_nao), nrow=m_nao, ncol=n_mc)

  bve_eq = c(bias2 = sum(rowMeans(w_est_eq - w_true)^2),
               var = sum(rowMeans((w_est_eq - rowMeans(w_est_eq))^2)),
               err = sum(rowMeans((w_est_eq - w_true)^2)))
  bve_uneq = c(bias2 = sum(rowMeans(w_est_uneq - w_true)^2),
               var = sum(rowMeans((w_est_uneq - rowMeans(w_est_uneq))^2)),
               err = sum(rowMeans((w_est_uneq - w_true)^2)))

  lst_nao[[paste(n_)]] = rbind(eq=bve_eq, uneq=bve_uneq)

  # ENSO
  w = replicate(n_mc, {
    data = sim(mu_enso, Sigma_enso, n_)
    w_eq = est_weights(data, 'equal')
    w_uneq = est_weights(data, 'unequal')
    c(w_eq, w_uneq)
  })

  w_est_eq = matrix(w[1, ], nrow=m_enso, ncol=n_mc, byrow=TRUE)
  w_est_uneq = w[-1, ]
  w_true = matrix(true_weights(Sigma_enso), nrow=m_enso, ncol=n_mc)

  bve_eq = c(bias2 = sum(rowMeans(w_est_eq - w_true)^2),
               var = sum(rowMeans((w_est_eq - rowMeans(w_est_eq))^2)),
               err = sum(rowMeans((w_est_eq - w_true)^2)))
  bve_uneq = c(bias2 = sum(rowMeans(w_est_uneq - w_true)^2),
               var = sum(rowMeans((w_est_uneq - rowMeans(w_est_uneq))^2)),
               err = sum(rowMeans((w_est_uneq - w_true)^2)))

  lst_enso[[paste(n_)]] = rbind(eq=bve_eq, uneq=bve_uneq)
}
bias_var_err_nao = lst_nao
bias_var_err_enso = lst_enso

save('bias-variance-weights.Rdata', list=c('bias_var_err_nao', 'bias_var_err_enso'))


