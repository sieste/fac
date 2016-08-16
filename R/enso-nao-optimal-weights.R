load('../data/enso-nao.Rdata')

# using full covariance
# ENSO
n_enso = nrow(enso)
Sigma_enso = var(enso)*(n_enso-1)/n_enso
mu_enso = colMeans(enso)

a_enso = mu_enso[1] - Sigma_enso[1, -1, drop=FALSE] %*% solve(Sigma_enso[-1, -1]) %*% mu_enso[-1]
w_enso = drop(Sigma_enso[1, -1, drop=FALSE] %*% solve(Sigma_enso[-1, -1]))
names(w_enso) = paste('w_', names(w_enso), sep='')
v_enso = Sigma_enso[1,1] - drop(Sigma_enso[1,-1,drop=FALSE] %*% solve(Sigma_enso[-1, -1]) %*% Sigma_enso[-1, 1, drop=FALSE])

out_enso = c(a=a_enso, w_enso, v=v_enso)
print(out_enso, digits=2)


# NAO
n_nao = ncol(nao)
Sigma_nao = var(nao)*(n_nao-1)/n_nao
mu_nao = colMeans(nao)

a_nao = mu_nao[1] - Sigma_nao[1, -1, drop=FALSE] %*% solve(Sigma_nao[-1, -1]) %*% mu_nao[-1]
w_nao = drop(Sigma_nao[1, -1, drop=FALSE] %*% solve(Sigma_nao[-1, -1]))
names(w_nao) = paste('w_', names(w_nao), sep='')
v_nao = Sigma_nao[1,1] - drop(Sigma_nao[1,-1,drop=FALSE] %*% solve(Sigma_nao[-1, -1]) %*% Sigma_nao[-1, 1, drop=FALSE])

out_nao = c(a=a_nao, w_nao, v=v_nao)
print(out_nao, digits=2)


# scale free, using correlation matrix
# using full covariance
# ENSO
n_enso = nrow(enso)
Sigma_enso = cor(enso)
mu_enso = rep(0, ncol(enso))

a_enso = mu_enso[1] - Sigma_enso[1, -1, drop=FALSE] %*% solve(Sigma_enso[-1, -1]) %*% mu_enso[-1]
w_enso = drop(Sigma_enso[1, -1, drop=FALSE] %*% solve(Sigma_enso[-1, -1]))
names(w_enso) = paste('w_', names(w_enso), sep='')
v_enso = Sigma_enso[1,1] - drop(Sigma_enso[1,-1,drop=FALSE] %*% solve(Sigma_enso[-1, -1]) %*% Sigma_enso[-1, 1, drop=FALSE])

out_enso = c(a=a_enso, w_enso, v=v_enso)
print(out_enso, digits=2)


# NAO
n_nao = nrow(nao)
Sigma_nao = cor(nao)
mu_nao = rep(0, ncol(nao))

a_nao = mu_nao[1] - Sigma_nao[1, -1, drop=FALSE] %*% solve(Sigma_nao[-1, -1]) %*% mu_nao[-1]
w_nao = drop(Sigma_nao[1, -1, drop=FALSE] %*% solve(Sigma_nao[-1, -1]))
names(w_nao) = paste('w_', names(w_nao), sep='')
v_nao = Sigma_nao[1,1] - drop(Sigma_nao[1,-1,drop=FALSE] %*% solve(Sigma_nao[-1, -1]) %*% Sigma_nao[-1, 1, drop=FALSE])

out_nao = c(a=a_nao, w_nao, v=v_nao)
print(out_nao, digits=2)
