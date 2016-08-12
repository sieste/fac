print('DOES NOT WORK PROPERLY AT THE MOMENT')

# I want an exponential distribution prior on the variance of the elements of
# the vector x - How should the elements of x be sampled?

p = 5 # length of x
C = 10 # parameter of the exponential distribution p(v) = C exp(-C*v)

sampls = matrix(nrow=0, ncol=p)

for (s in 1:10000) {
  x = vector()
  m = 0
  for (i in 1:p) {
    tmp = rnorm(1, m, p/sqrt(2*C*(p-1)))
    x = c(x, tmp)
    m = m + tmp / (p - 1)
  }
  sampls = rbind(sampls, x)
}

v = apply(sampls, 1, function(x) mean((x-mean(x))^2))
plot(density(v, from=0, to=10/C))
xx = seq(0, 10/C, length.out=100)
lines(xx, C*exp(-xx*C), lty=2)

