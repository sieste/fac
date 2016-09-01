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
  xx_uneq = nvec+1
  yy_eq = rbind(data[, 'mean_eq'],
                data[, 'mean_eq'] - data[, 'sd_eq'],
                data[, 'mean_eq'] + data[, 'sd_eq'])
  yy_uneq = rbind(data[, 'mean_uneq'],
                  data[, 'mean_uneq'] - data[, 'sd_uneq'],
                  data[, 'mean_uneq'] + data[, 'sd_uneq'])
  
  plot(NULL, xlim=range(c(xx_eq, xx_uneq)), ylim=range(c(yy_eq, yy_uneq)), 
       axes=FALSE, ann=FALSE)
  lines(xx_eq, yy_eq[1, ])
  segments(xx_eq, yy_eq[2, ], xx_eq, yy_eq[3, ])
  lines(xx_uneq, yy_uneq[1, ], lty=2)
  segments(xx_uneq, yy_uneq[2, ], xx_uneq, yy_uneq[3, ], lty=2)

  legend('topright', c('equal weights', 'unequal weights'), lty=c(1,2))
  axis(1)
  axis(2, las=2)
  mtext(side=1, text='number of training data', line=3)
  mtext(side=2, text='score (mean +/- stdev)', line=3)
  title(main=label)
  box()

  
}





plot(NULL, xlim=c(0,200), ylim=c(.3,.7),  ann=FALSE, axes=FALSE)
abline(h=.5, col='gray')
lines(nvec, impr_nao[,'se'])
lines(nvec, impr_nao[,'crps'], lty=2)
lines(nvec, impr_nao[,'logs'], lty=3)
legend('topleft', c('sq. error', 'CRPS', 'LOGS'), lty=1:3)
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
legend('topleft', c('sq. error', 'CRPS', 'LOGS'), lty=1:3)
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
