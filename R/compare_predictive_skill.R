source('functions.R')

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

print('one factor:')
print(fa_loo(x, y, nf=1, Psi.lower=c(1e-6, Psi.lower))) 
print('two factors:')
print(fa_loo(x, y, nf=2, Psi.lower=c(1e-6, Psi.lower))) 

print('one factor spherical:')
print(fa_loo(x, y, nf=1, sph=TRUE))
print('two factors: spherical:')
print(fa_loo(x, y, nf=2, sph=TRUE))

# linear regression of the mmem
n = length(y)
m = rowMeans(x)
score = 0
for (i in 1:n) {
  
  mm = m[-i]
  yy = y[-i]

  M_p = mean(yy) + cov(mm,yy) / var(mm) * (m[i] - mean(mm))
  C_p = var(yy) - cov(mm,yy)^2 / var(mm)

  # increment leave-one-out score
  score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)
}
print('simple linear regression:')
print(score/n)


# multiple linear regression
n = length(y)
X = cbind(y, x)
score = 0
for (i in 1:n) {
  
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
print('multiple linear regression:')
print(score/n)


# ridge multiple regression
n = length(y)
score = 0
for (i in 1:n) {
  
  xx = x[-i, ]
  yy = y[-i]
  mod = cv.glmnet(xx, yy, alpha=0, grouped=FALSE)
  beta = as.vector(coef(mod, s=mod$lambda.min))
  C_p = min(mod$cvm)
  M_p = predict(mod, s=mod$lambda.min, x[i,,drop=FALSE])
  
  # increment leave-one-out score
  score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)

}
print('multiple ridge regression:')
print(score/n)



# ridge lasso regression
n = length(y)
score = 0
for (i in 1:n) {
  
  xx = x[-i, ]
  yy = y[-i]
  mod = cv.glmnet(xx, yy, alpha=1, grouped=FALSE)
  beta = as.vector(coef(mod, s=mod$lambda.min))
  C_p = min(mod$cvm)
  M_p = predict(mod, s=mod$lambda.min, newx=x[i,,drop=FALSE])
  
  # increment leave-one-out score
  score = score - dnorm(y[i], M_p, sqrt(C_p), log=TRUE)

}
print('multiple lasso regression:')
print(score/n)


# principal component regression
library(pls)
score1 = 0
score2 = 0
X = data.frame(cbind(obs=y, x))
for (i in 1:n) {
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
print('principal component regression:')
print(score1/n)
print(score2/n)



