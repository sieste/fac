source('functions.R')

# ###############################################################
# # EL NINO DATA from CAIO
# ###############################################################
# 
# # factor analysis of ensemble means
# nino = list(
#   cfs = read.table('../data/cfsv2_tos_nino34_data.txt') - 273.15,
#   cmc = read.table('../data/cmc2_tos_nino34_data.txt') - 273.15,
#   gfdl = read.table('../data/gfdl_tos_nino34_data.txt'),
#   mf = read.table('../data/mf3_tos_nino34_data.txt') - 273.15,
#   nasa = read.table('../data/nasa_tos_nino34_data.txt'),
#   ec = read.table('../data/syst4_tos_nino34_data.txt') - 273.15
# )
# 
# x = sapply(nino, rowMeans)
# y = drop(as.matrix(read.table('../data/obs_tos_nino34_data.txt') - 273.15))
# Psi.lower = sapply(nino, function(x) mean(apply(x, 1, var)/ncol(x)))


# ###############################################################
# # DEMETER NAO data
# ###############################################################
#
# nao = list(
# #  cerfacs = read.table('../data/demeter_nao_cerfacs.dat'),
# #  ingv = read.table('../data/demeter_nao_ingv.dat'),
#   ecmwf = read.table('../data/demeter_nao_ecmwf.dat'),
#   lodyn = read.table('../data/demeter_nao_lodyn.dat'),
#   metfr = read.table('../data/demeter_nao_metfr.dat'),
#   mpi = read.table('../data/demeter_nao_mpi.dat'),
#   ukmo = read.table('../data/demeter_nao_ukmo.dat'),
#   obs = read.table('../data/obs_nao_climexp.txt')
# )
# 
# nao_df = data.frame(year=1995)
# for (i in 1:length(nao)) {
#   df_tmp = data.frame(nao[[i]][,1], rowMeans(nao[[i]][, -1, drop=FALSE]))
#   names(df_tmp) = c('year', names(nao)[i])
#   nao_df = merge(nao_df, df_tmp, by='year', all=TRUE)
# }
# 
# nao_df = nao_df[!is.na(rowSums(nao_df)), ]
# x = as.matrix(nao_df[-which(names(nao_df) %in% c('year','obs'))])
# y = drop(as.matrix(nao_df['obs']))
# Psi.lower = rep(1e-6, ncol(x))

###############################################################
# CAIO's NAO data
###############################################################

nao = list(
  cm1 = cbind(1982:2010, read.table('~/folders/nino34-combination/nao-test/caio-data/CMC2-CanCM4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  cm2 = cbind(1982:2010, read.table('~/folders/nino34-combination/nao-test/caio-data/CMC1-CanCM3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  ecmwf = cbind(1982:2010, read.table('~/folders/nino34-combination/nao-test/caio-data/ECMWFS4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
  metfr = cbind(1982:2010, read.table('~/folders/nino34-combination/nao-test/caio-data/MFS3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')),
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


Cvec = seq(2, 6, .2)
score = c()
for (C in Cvec) {
  score = c(score, mean(fa_loo(x, y, fa_ple, nf=1L, C=10^C)))
}

plot(Cvec, score)
