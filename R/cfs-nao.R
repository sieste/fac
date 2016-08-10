load('../data/cfs-psl.Rdata')
print(call)
# "cfs = loadECOMS(dataset='CFSv2_seasonal', var='psl', lonLim=c(-90,60),
# latLim=c(20,90), season=c(12,1,2), leadMonth=1, time='DD', aggr.d='mean',
# aggr.m='mean')"

nao = cfs$Data

# average over longitude (4th dimension)
nao = apply(X=nao, MARGIN=1:3, FUN=mean)

# weighted averages over latitude (3rd dimension)
w.north = cos((55:90) / 360 * 2 * pi)
north   = apply(X=nao[,,36:71], MARGIN=c(1,2), FUN=weighted.mean, w=w.north)
w.south = cos((20:55) / 360 * 2 * pi)
south   = apply(X=nao[,,1:36], MARGIN=c(1,2), FUN=weighted.mean, w=w.south)
nao     = south - north

# aggregate months into seasons
yy = rep(1:27, each=3)
nao = apply(nao, 1, function(member) tapply(member, yy, mean))

# put NAs for missing years 1982 and 2010
nao = rbind(NA_real_, nao, NA_real_)

# save
write.table(nao,"../data/CFS2-NAO-hind-ens-members-ic-Nov-val-DJF-1982-2010.txt",row.names=FALSE,col.names=FALSE)

