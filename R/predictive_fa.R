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

print(fa_loo(x, y, nf=1, Psi.lower=c(1e-6, Psi.lower))) 
print(fa_loo(x, y, nf=2, Psi.lower=c(1e-6, Psi.lower))) 
print(fa_loo(x, y, nf=3, Psi.lower=c(1e-6, Psi.lower))) 
print(fa_loo(x, y, nf=4, Psi.lower=c(1e-6, Psi.lower))) 



print('')

###############################################################
# DEMETER NAO data
###############################################################

nao = list(
  cerfacs = read.table('../data/demeter_nao_cerfacs.dat'),
  ingv = read.table('../data/demeter_nao_ingv.dat'),
  mpi = read.table('../data/demeter_nao_mpi.dat'),
  lodyn = read.table('../data/demeter_nao_lodyn.dat'),
  ecmwf = read.table('../data/demeter_nao_ecmwf.dat'),
  metfr = read.table('../data/demeter_nao_metfr.dat'),
  ukmo = read.table('../data/demeter_nao_ukmo.dat'),
  obs = read.table('../data/obs_nao_climexp.txt')
)

nao_df = data.frame(year=1995)
for (i in 1:length(nao)) {
  df_tmp = data.frame(nao[[i]][,1], rowMeans(nao[[i]][, -1, drop=FALSE]))
  names(df_tmp) = c('year', names(nao)[i])
  nao_df = merge(nao_df, df_tmp, by='year', all=TRUE)
}

nao_df = nao_df[!is.na(rowSums(nao_df)), -1]
y = c(as.matrix(nao_df['obs']))
x = as.matrix(nao_df[names(nao_df) != 'obs'])

print(fa_loo(x, y, 1, Psi.lower=1e-1))
print(fa_loo(x, y, 2, Psi.lower=1e-1))
print(fa_loo(x, y, 3, Psi.lower=1e-1))
print(fa_loo(x, y, 4, Psi.lower=1e-1))

# optimum is for one factor for most combination, but two factors for x=(ecmwf,
# metfr, ukmo)



