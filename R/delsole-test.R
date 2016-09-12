load('../data/enso-nao.Rdata')
uneq_test = function(y, x) {
  R2_eq = summary(lm(y~rowMeans(x)))$r.squared
  R2_uneq = summary(lm(y~x))$r.squared
  N = length(y)
  M = ncol(x)
  delta = (R2_uneq - R2_eq) / (1 - R2_eq)
  Fstat = delta / (1 - delta) * (N - M - 1) / (M - 1)
  pval = pf(Fstat, N-1, N-M-1, lower.tail=FALSE)
  return(c(Fstat, pval))
}

uneq_test(y=enso[,1], x=enso[,-1])
uneq_test(y=nao[,1], x=nao[,-1])

# in both data sets, the test does not lead to rejection of the null hypothesis
# of equal weights




# monte carlo experiment to check how well the delsole test predicts whether
# equal or unequal weighting is better
delsole_test_monte_carlo_loop = 
function(Mu, S, n_mc = 1e4, 
         n_vec = c(10, 20, 50, 100)) 
{

  p = length(Mu)
  res = array(NA_real_, 
              dim=c(length(n_vec), n_mc, 3), 
              dimnames=list(
                  paste(n_vec), 
                  NULL, 
                  c('se_eq', 'se_uneq', 'pval_delsole')
              ))


  for (i_n in seq_along(n_vec)) {
    for (i_mc in seq_len(n_mc)) {
      n_ = n_vec[i_n]

      if (i_mc %% 1e4 == 0) {
        print(match.call())
        cat('n = ', n_, ' i_mc = ', i_mc, '\n')
      }
  
      # simulate independent training data of size n, and one test sample
      data = rmvnorm(n_+1, mean=Mu, sigma=S)
      y_tr = data[(1:n_), 1]
      x_tr = cbind(1, data[(1:n_), -1])
      x2_tr = cbind(1, rowMeans(data[(1:n_), -1]))
    
      # fit unequal weights
      XtXm1_uneq = solve(crossprod(x_tr))
      beta_uneq = XtXm1_uneq %*% crossprod(x_tr, y_tr)
      rss_uneq = sum((y_tr - x_tr %*% beta_uneq)^2)
    
      # fit equal weights
      XtXm1_eq = solve(crossprod(x2_tr))
      beta_eq = XtXm1_eq %*% crossprod(x2_tr, y_tr)
      rss_eq = sum((y_tr - x2_tr %*% beta_eq)^2)
    
      # squared prediction errors
      y_te = data[n_ + 1, 1]
      x_te = c(1, data[n_ + 1, -1])
      x2_te = c(1, mean(data[n_ + 1, -1]))
      se_uneq = (y_te - sum(x_te * beta_uneq))^2
      se_eq = (y_te - sum(x2_te * beta_eq))^2
    
      # DelSole's F-test: eq 16 of DelSole (2013) "Unequal and Equal Weighting"
      M_ = length(Mu) - 1
      Fstat = (rss_eq - rss_uneq) / rss_uneq * (n_ - M_ - 1) / (M_ - 1)
      pval = pf(Fstat, n_-1, n_-M_-1, lower.tail=FALSE)
       
      res[i_n, i_mc, ] = c(se_eq, se_uneq, pval)
    }
  }

  return(res)

}

# some random number experiments with the enso data
n_enso = nrow(enso)
S = var(enso) * (n_enso-1) / n_enso
Mu = colMeans(enso)
res_enso = delsole_test_monte_carlo_loop(Mu=Mu, S=S)

n_nao = nrow(nao)
S_nao = var(nao) * (n_nao-1) / n_nao
mu_nao = colMeans(nao)
res_nao = delsole_test_monte_carlo_loop(mu_nao, S_nao)

# calculate mean squared out-of-sample errors

# check how often equal weighting is worse than unequal weighting

# power: check how often the delsole test yields a pvalue less than 0.05


