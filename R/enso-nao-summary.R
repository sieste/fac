# ENSO data
enso = list(
  obs = as.matrix(read.table('../data/obs_tos_nino34_data.txt') - 273.15),
  cfs = read.table('../data/cfsv2_tos_nino34_data.txt') - 273.15,
  cmc = read.table('../data/cmc2_tos_nino34_data.txt') - 273.15,
  gfdl = read.table('../data/gfdl_tos_nino34_data.txt'),
  mf = read.table('../data/mf3_tos_nino34_data.txt') - 273.15,
  nasa = read.table('../data/nasa_tos_nino34_data.txt'),
  ec = read.table('../data/syst4_tos_nino34_data.txt') - 273.15
)

enso = sapply(enso, rowMeans)
n = nrow(enso)

cat('ENSO MEANS:\n')
print(colMeans(enso), digits=4)
cat('ENSO STDEVS:\n')
print(apply(enso, 2, sd)*(n-1)/n, digits=2)
cat('ENSO CORRELATIONS:\n')
print(cor(enso), digits=2)
cat('\n')


# DEMETER NAO data
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
nao = cbind(obs=y, x)


cat('NAO MEANS:\n')
print(colMeans(nao), digits=2)
cat('NAO STDEVS:\n')
print(apply(nao, 2, sd)*(n-1)/n, digits=2)
cat('NAO CORRELATIONS:\n')
print(cor(nao), digits=1)
cat('\n')


# save for later use (e.g. to simulate similar data)
save(file='../data/enso-nao.Rdata', list=c('enso', 'nao'))

