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
# scale models by 1/1000
nao[, -1] = nao[, -1] / 1000
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

# equal weighting
M_enso = ncol(enso)-1
A = rbind(c(1, rep(0, M_enso)), c(0, rep(1/M_enso, M_enso)))
Sigma_enso_mmm = A %*% Sigma_enso %*% t(A)
mu_enso_mmm = A %*% mu_enso

a_enso_mmm = mu_enso_mmm[1] - Sigma_enso_mmm[1,2] / Sigma_enso_mmm[2,2] * mu_enso_mmm[2]
w_enso_mmm = Sigma_enso_mmm[1,2] / Sigma_enso_mmm[2,2]
s_enso_mmm = Sigma_enso_mmm[1,1] - Sigma_enso_mmm[1,2]^2 / Sigma_enso_mmm[2,2]

cat('enso equal weighting: a, w, s:\n')
cat(a_enso_mmm, w_enso_mmm, s_enso_mmm, '\n')


# equal weighting
M_nao = ncol(nao)-1
A = rbind(c(1, rep(0, M_nao)), c(0, rep(1/M_nao, M_nao)))
Sigma_nao_mmm = A %*% Sigma_nao %*% t(A)
mu_nao_mmm = A %*% mu_nao

a_nao_mmm = mu_nao_mmm[1] - Sigma_nao_mmm[1,2] / Sigma_nao_mmm[2,2] * mu_nao_mmm[2]
w_nao_mmm = Sigma_nao_mmm[1,2] / Sigma_nao_mmm[2,2]
s_nao_mmm = Sigma_nao_mmm[1,1] - Sigma_nao_mmm[1,2]^2 / Sigma_nao_mmm[2,2]

cat('nao equal weighting: a, w, s:\n')
cat(a_nao_mmm, w_nao_mmm, s_nao_mmm, '\n')
