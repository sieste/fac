# download from climate explorer
nao = read.table('http://climexp.knmi.nl/data/inao.dat', 
                 na.strings='-999.9000', row.names=1)
n = nrow(nao)

# average DJF nao
djf = rowMeans(cbind(nao[1:(n-1), 12], nao[2:n, 1], nao[2:n, 2]))
names(djf) = rownames(nao)[2:n]
write.table(round(djf, 2), '../data/obs_nao_climexp.txt', col.names=FALSE)


