source('functions.R')



compare_skill = function(x,y) {


  #################################################################
  print('factor analysis, one factor:')
  print(fa_loo(x, y, fa_mle, nf=1L, Psi.lower=c(1e-6, Psi.lower))) 

  print('factor analysis, two factors:')
  print(fa_loo(x, y, fa_mle, nf=2L, Psi.lower=c(1e-6, Psi.lower))) 
  

  #################################################################
  print('factor analysis, one factor, spherical noise:')
  print(fa_loo(x, y, fa_mle_sph, nf=1L))

  print('factor analysis, two factors, spherical noise:')
  print(fa_loo(x, y, fa_mle_sph, nf=2L))
  

  #################################################################
  # linear regression of the mmem
  n = length(y)
  m = rowMeans(x)
  score = 0
  for (i in 1L:n) {
    
    mm = m[-i]
    yy = y[-i]
  
    M_p = mean(yy) + cov(mm,yy) / var(mm) * (m[i] - mean(mm))
    C_p = var(yy) - cov(mm,yy)^2 / var(mm)
  
    # increment leave-one-out score
    score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)
  }
  print('simple linear regression on the mme mean:')
  print(score/n)
  
  
  #################################################################
  # multiple linear regression
  n = length(y)
  X = cbind(y, x)
  score = 0
  for (i in 1L:n) {
    
      # mle of covariance matrix 
      M = colMeans(X[-i, ])
      C = var(X[-i, ])
  
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
  print('multiple linear regression on the individual ensemble means:')
  print(score/n)
  
  
  #################################################################
  # ridge multiple regression
  library(glmnet)
  n = length(y)
  score = 0
  for (i in 1L:n) {
    
    xx = x[-i, ]
    yy = y[-i]
    mod = cv.glmnet(xx, yy, alpha=0, grouped=FALSE)
    beta = as.vector(coef(mod, s=mod$lambda.min))
    C_p = min(mod$cvm)
    M_p = predict(mod, s=mod$lambda.min, x[i,,drop=FALSE])
    
    # increment leave-one-out score
    score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)
  
  }
  print('penalised multiple regression, ridge penalty:')
  print(score/n)
  
  
  
  #################################################################
  # ridge lasso regression
  n = length(y)
  score = 0
  for (i in 1L:n) {
    
    xx = x[-i, ]
    yy = y[-i]
    mod = cv.glmnet(xx, yy, alpha=1, grouped=FALSE)
    beta = as.vector(coef(mod, s=mod$lambda.min))
    C_p = min(mod$cvm)
    M_p = predict(mod, s=mod$lambda.min, newx=x[i,,drop=FALSE])
    
    # increment leave-one-out score
    score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)
  
  }
  print('penalised multiple regression, lasso penality:')
  print(score/n)
  
  #################################################################
  # principal component regression
  library(pls)
  score1 = 0
  score2 = 0
  X = data.frame(cbind(obs=y, x))
  for (i in 1L:n) {
    xx = X[-i, ]
    mod = pcr(formula=obs~., ncomp=2, data=xx)
  
    C_p1 = mean((y[-i] - predict(mod, comps=1))^2)
    M_p1 = predict(mod, newdata=X[i,], comps=1)
    C_p2 = mean((y[-i] - predict(mod, comps=1:2))^2)
    M_p2 = predict(mod, newdata=X[i,], comps=1:2)
    
    # increment leave-one-out score
    score1 = score1 - dnorm(y[i], M_p1, sqrt(C_p1), log=TRUE)
    score2 = score2 - dnorm(y[i], M_p2, sqrt(C_p2), log=TRUE)
  
  }
  print('principal component regression, one component:')
  print(score1/n)
  print('principal component regression, two components:')
  print(score2/n)


}




###############################################################
# EL NINO DATA from CAIO
###############################################################

# factor analysis of ensemble means
nino = list(
  cfs = read.table('../data/cfsv2_tos_nino34_data.txt') - 273.15,
  cmc = read.table('../data/cmc2_tos_nino34_data.txt') - 273.15,
  gfdl = read.table('../data/gfdl_tos_nino34_data.txt'),
  mf = read.table('../data/mf3_tos_nino34_data.txt') - 273.15,
  nasa = read.table('../data/nasa_tos_nino34_data.txt'),
  ec = read.table('../data/syst4_tos_nino34_data.txt') - 273.15
)

x = sapply(nino, rowMeans)
y = drop(as.matrix(read.table('../data/obs_tos_nino34_data.txt') - 273.15))
Psi.lower = sapply(nino, function(x) mean(apply(x, 1, var)/ncol(x)))

print('################################################')
print('seasonal El Nino forecasts')
compare_skill(x,y)


###############################################################
# DEMETER NAO data
###############################################################

nao = list(
#  cerfacs = read.table('../data/demeter_nao_cerfacs.dat'),
#  ingv = read.table('../data/demeter_nao_ingv.dat'),
  ecmwf = read.table('../data/demeter_nao_ecmwf.dat'),
  lodyn = read.table('../data/demeter_nao_lodyn.dat'),
  metfr = read.table('../data/demeter_nao_metfr.dat'),
  mpi = read.table('../data/demeter_nao_mpi.dat'),
  ukmo = read.table('../data/demeter_nao_ukmo.dat'),
  obs = read.table('../data/obs_nao_climexp.txt')
)

nao_df = data.frame(year=1995)
for (i in 1:length(nao)) {
  df_tmp = data.frame(nao[[i]][,1], rowMeans(nao[[i]][, -1, drop=FALSE]))
  names(df_tmp) = c('year', names(nao)[i])
  nao_df = merge(nao_df, df_tmp, by='year', all=TRUE)
}

nao_df = nao_df[!is.na(rowSums(nao_df)), ]
x = as.matrix(nao_df[-which(names(nao_df) %in% c('year','obs'))])
y = drop(as.matrix(nao_df['obs']))
Psi.lower = rep(1e-6, ncol(x))

print('################################################')
print('Demeter seasonal NAO forecasts')
compare_skill(x,y)



###############################################################
# CAIO's seasonal NAO data (1982 - 2010)
###############################################################

nao = list(
  cm1 = cbind(1982:2010, read.table('../data/CMC2-CanCM4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  cm2 = cbind(1982:2010, read.table('../data/CMC1-CanCM3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  ecmwf = cbind(1982:2010, read.table('../data/ECMWFS4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  metfr = cbind(1982:2010, read.table('../data/MFS3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  obs = read.table('../data/obs_nao_climexp.txt')
)

nao_df = data.frame(year=1995)
for (i in 1:length(nao)) {
  df_tmp = data.frame(nao[[i]][,1], rowMeans(nao[[i]][, -1, drop=FALSE]))
  names(df_tmp) = c('year', names(nao)[i])
  nao_df = merge(nao_df, df_tmp, by='year', all=TRUE)
}
nao_df = nao_df[!is.na(rowSums(nao_df)), ]
x = as.matrix(nao_df[-which(names(nao_df) %in% c('year','obs'))])
y = drop(as.matrix(nao_df['obs']))
Psi.lower = rep(1e-6, ncol(x))

print('################################################')
print('more recent seasonal NAO forecasts')
compare_skill(x,y)
