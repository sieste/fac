load('n-dependence-scores-summary.Rdata')

pdf('n-dependence.pdf')

lst = list(
NAO_log_score=logs_nao,
NAO_crps=crps_nao,
NAO_squared_error=se_nao,
ENSO_log_score=logs_enso,
ENSO_crps=crps_enso,
ENSO_squared_error=se_enso)

for (i in seq_along(lst)) {
  data = lst[[i]]
  label = names(lst)[i]
  xx_eq = nvec
  xx_uneq = nvec
  xx_clim = nvec
  xx_fa = nvec
  mean_eq = data[, 'mean_eq']
  mean_uneq = data[, 'mean_uneq']
  mean_clim = data[, 'mean_cl']
  mean_fa = data[, 'mean_fa']
  sd_eq = data[, 'sd_eq']
  sd_uneq = data[, 'sd_uneq']
  sd_clim = data[, 'sd_cl']
  sd_fa = data[, 'sd_fa']


  plot(NULL, xlim=range(nvec), 
             ylim=range(c(mean_eq, mean_uneq, mean_clim, mean_fa)), 
             axes=FALSE, ann=FALSE)
  lines(xx_clim, mean_clim, lty=3)
  lines(xx_eq, mean_eq)
  lines(xx_uneq, mean_uneq, lty=2)
  lines(xx_fa, mean_fa, lty=4)

  legend('topright', c('equal weights', 'unequal weights', 'climatology', 'factor analysis'), lty=c(1,2,3,4))
  axis(1)
  axis(2, las=2)
  mtext(side=1, text='number of years for training', line=3)
  mtext(side=2, text='mean score', line=3)
  title(main=paste(label, 'mean'))
  box()

  plot(NULL, xlim=range(nvec), 
             ylim=range(c(sd_eq, sd_uneq, sd_clim, sd_fa)), 
             axes=FALSE, ann=FALSE)
  lines(xx_clim, sd_clim, lty=3)
  lines(xx_eq, sd_eq)
  lines(xx_uneq, sd_uneq, lty=2)
  lines(xx_fa, sd_fa, lty=4)

  legend('topright', c('equal weights', 'unequal weights', 'climatology', 'factor analysis'), lty=c(1,2,3,4))
  axis(1)
  axis(2, las=2)
  mtext(side=1, text='number of years for training', line=3)
  mtext(side=2, text='sd score', line=3)
  title(main=paste(label,'stdev'))
  box()

  
}





plot(NULL, xlim=c(0,200), ylim=c(.3,.7),  ann=FALSE, axes=FALSE)
abline(h=.5, col='gray')
lines(nvec, impr_nao[,'se'])
lines(nvec, impr_nao[,'crps'], lty=2)
lines(nvec, impr_nao[,'logs'], lty=3)
legend('topleft', c('SQERR', 'CRPS', 'LOGS'), lty=1:3)
axis(1)
axis(2, las=2)
mtext(side=1, text='number of training data', line=3)
mtext(side=2, 'probability that unequal weighting improves the forecast', line=3)
title(main='NAO: P(score_unequal < score_equal)')
box()

plot(NULL, xlim=c(0,200), ylim=c(.3,.7), ann=FALSE, axes=FALSE)
abline(h=.5, col='gray')
lines(nvec, impr_enso[,'se'])
lines(nvec, impr_enso[,'crps'], lty=2)
lines(nvec, impr_enso[,'logs'], lty=3)
legend('topleft', c('SQERR', 'CRPS', 'LOGS'), lty=1:3)
axis(1)
axis(2, las=2)
mtext(side=1, text='number of training data', line=3)
mtext(side=2, 'probability that unequal weighting improves the forecast', line=3)
title(main='ENSO: P(score_unequal < score_equal)')
box()

# illustration
set.seed(123)
n = 20
xx = 1:n
s = rnorm(n)
obs = s + rnorm(n, 0, .5) + .4
ens = s + matrix(rnorm(n*4, 0, .5), n, 4)
plot(NULL, xlim=range(1, n+5), ylim=range(c(ens,obs)), ann=FALSE, axes=FALSE)
lines(xx, obs, col=gray(.5), lwd=7)
matlines(xx, ens, col='black', lty=1)

abline(v=n+1, lty=2)

ens2 = rnorm(1) + rnorm(4, 0, .5)
matlines(c(n, n+2), rbind(ens[n,], ens2), col='black', lty=1)
points(rep(n+2, 4), ens2)
lines(n+2+2*dnorm(seq(-1,1.8,.1), mean=.4, sd=.4), seq(-1,1.8,.1))
points(n+2, rnorm(1)+1, pch=15, cex=2, col=gray(.5))


dev.off()

