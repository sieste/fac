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
with(fa_mle(x, 1), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 2), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 3), print(c(AIC=AIC, BIC=BIC)))

###############################################################
# AIC prefers nf=2, BIC prefers nf=1 (same if obs is included)
###############################################################



###############################################################
# NAO DATA from CAIO
###############################################################

# factor analysis of ensemble means
nao = list(
  cfs = read.table('../data/CFS2-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt'),
  cmc1 = read.table('../data/CMC1-CanCM3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt'), 
  cmc2 = read.table('../data/CMC2-CanCM4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt'),
  ec = read.table('../data/ECMWFS4-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt'),
  mf = read.table('../data/MFS3-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt')
)

x = sapply(nao, rowMeans)
x = x[-c(1,29), ]

with(fa_mle(x, 1), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 2), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 3), print(c(AIC=AIC, BIC=BIC)))

###############################################################
# AIC prefers nf=2, BIC prefers nf=1 
###############################################################




###############################################################
# DEMETER NAO data
###############################################################

nao = list(
  cerfacs = read.table('../data/demeter_nao_cerfacs.dat'),
  ingv = read.table('../data/demeter_nao_ingv.dat'),
  ecmwf = read.table('../data/demeter_nao_ecmwf.dat'),
  lodyn = read.table('../data/demeter_nao_lodyn.dat'),
  metfr = read.table('../data/demeter_nao_metfr.dat'),
  mpi = read.table('../data/demeter_nao_mpi.dat'),
  ukmo = read.table('../data/demeter_nao_ukmo.dat')
)

nao_df = data.frame(year=1995)
for (i in 1:length(nao)) {
  df_tmp = data.frame(nao[[i]][,1], rowMeans(nao[[i]][, -1]))
  names(df_tmp) = c('year', names(nao)[i])
  nao_df = merge(nao_df, df_tmp, by='year', all=TRUE)
}

x = nao_df[!is.na(rowSums(nao_df)), -1]

source('functions.R')
with(fa_mle(x, 1), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 2), print(c(AIC=AIC, BIC=BIC)))
with(fa_mle(x, 3), print(c(AIC=AIC, BIC=BIC)))

###############################################################
# the 1 factor model is favored by both AIC and BIC, and for various
# combinations of models
###############################################################


