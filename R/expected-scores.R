expected_scores = function(Sigma) {
  Sigma = as.matrix(Sigma)
  if (ncol(Sigma)==1) {
    s_yx = drop(Sigma)
  } else {
    s_yx = drop(Sigma[1,1] - Sigma[1, -1, drop=FALSE] %*% solve(Sigma[-1, -1]) %*% Sigma[-1, 1, drop=FALSE])
  }
  return(c(mse=s_yx, crps=sqrt(s_yx/pi), ign=.5*(log(2*pi)+log(s_yx)+1)))
}

load('../data/enso-nao.Rdata')
n_enso = nrow(enso)
p_enso = ncol(enso)
M_enso = p_enso - 1
S_enso = var(enso) * (n_enso-1) / n_enso
A_enso = rbind(c(1,rep(0, M_enso)), c(0, rep(1/M_enso, M_enso)))

S_enso_clim = S_enso[1,1]
S_enso_mmm = A_enso %*% S_enso %*% t(A_enso)

scores_enso = rbind(
  clim= expected_scores(S_enso_clim),
  equal = expected_scores(S_enso_mmm),
  unequal = expected_scores(S_enso)
  )

cat('ENSO\n')
print(scores_enso, digits=2)
cat('\n')

n_nao = nrow(nao)
p_nao = ncol(nao)
M_nao = p_nao - 1
S_nao = var(nao) * (n_nao-1) / n_nao
A_nao = rbind(c(1,rep(0, M_nao)), c(0, rep(1/M_nao, M_nao)))

S_nao_clim = S_nao[1,1]
S_nao_mmm = A_nao %*% S_nao %*% t(A_nao)

scores_nao = rbind(
  clim= expected_scores(S_nao_clim),
  equal = expected_scores(S_nao_mmm),
  unequal = expected_scores(S_nao)
  )

cat('NAO\n')
print(scores_nao, digits=5)
cat('\n')


cat('\n')
cat('skill scores:\n\n')
cat('ENSO\n')
print((scores_enso[1,]-scores_enso[2,])/(scores_enso[1,]-scores_enso[3,]), digits=2)
cat('\n')
cat('NAO\n')
print((scores_nao[1,]-scores_nao[2,])/(scores_nao[1,]-scores_nao[3,]), digits=2)

