# Selection

# Change in allele frequency
# changes in allele frequency
t <- 1:500000
f0 <- 0.01
s <- 0.00001
plot(f0/(f0+(1-s)^t*(1-f0)), ylab="frequency", xlab="generations")

# exponential distribution
lines(f0/(f0+exp(-s*t)*(1-f0)), type="l", col="red", lwd=2)

# Special cases
## directional selection
s <- 0.1 # selection coefficient
## additive
f <- rep(0,1000)
f[1] <- 0.01
for (t in 2:1000) f[t] <- f[t-1] + s*f[t-1]*(1-f[t-1])
plot(f, type="l", col="red")
legend("bottomright", col=c("red","black","blue"), legend=c("additive",
                                  "dominant", "recessive"), lty=1, lwd=2)
## dominant
f <- rep(0,1000)
f[1] <- 0.01
for (t in 2:1000) f[t] <- f[t-1] + s*f[t-1]*(1-f[t-1])^2 / (1 - s*(1-f[t-1]^2))
lines(f, type="l", col="black", lwd=2)

## recessive
f <- rep(0,1000)
f[1] <- 0.01
for (t in 2:1000) f[t] = f[t-1] + (s*(f[t-1])^2*(1-f[t-1])) / (1 - s*(2*f[t-1]*(1-f[t-1]) + (1-f[t-1])^2))
lines(f, type="l", col="blue", lwd=2)

## Selection and drift
plotTrajectory(simulateTrajectory(s=0.001, N=100, t=100, nrepl=100))

## change "s" and "N"!
