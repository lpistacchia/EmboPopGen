#Expected allele frequency
#Genetic Drift

################
#1 Calculate allele frequency after one binomial sampling (over one generation)
N <- 50
fA <- 0.53
frequencies <- mean(rbinom(1, 2*N, fA) / (2*N))
print(frequencies)
################

#2 
N <- 150
fA <- 0.53
gen <- 10
nrepl <- 20

N <- 100
gen <- 1000 
nrepl <- 20 
for (j in 1:nrepl) {
  fA <- rep(NA, gen)
  fA[1] <- 0.5
  for (t in 1:(gen-1)) fA[t+1] <- rbinom(1, 2*N, fA[t]) / (2*N)
  if (j==1) plot(x=1:gen, y=fA, xlab="generations", type="l", ylim=c(0,1), col=rainbow(nrepl)[j]) else lines(x=1:gen, y=fA, col=rainbow(nrepl)[j])
}

################

#2 Initialize matrix to store frequencies for many generations
frequencies <- matrix(0, nrepl, gen)
for (j in 1:nrepl) {
  current_fA <- fA
  for (i in 1:gen) {
    # Simulate allele frequency for the next generation
    current_fA <- rbinom(1, 2*N, current_fA) / (2*N)
    frequencies[j, i] <- current_fA
  }
}

print(frequencies)
hist(frequencies)


################


#3
# small population (blue)
N <- 100
for (j in 1:10) {
  fA <- rep(NA, 100); fA[1] <- 0.5
  for (t in 1:99) fA[t+1] <- rbinom(1, 2*N, fA[t]) / (2*N)
  if (j==1) plot(x=1:100, y=fA, type="l", ylim=c(0,1), col="blue") else lines(x=1:100, y=fA, col="blue")
}

# large population (red)
N <- 1000
for (j in 1:10) {
  fA <- rep(NA, 100); fA[1] <- 0.5
  for (t in 1:99) fA[t+1] <- rbinom(1, 2*N, fA[t]) / (2*N)
  lines(x=1:100, y=fA, col="red")
}


rbinom(1, 20, 0.5) / 20
rbinom(1, 20000, 0.5) / 20000

# simulateTrajectory
simulateTrajectory <- function(s, N, t=500, nrepl=100) {
  
  cat("2Ns =",2*N*s,"\n")
  
  # initialise frequencies
  fA <- matrix(NA, nrow=nrepl, ncol=t)
  fA[,1] <- 1/(2*N)
  
  # viability
  vAA <- 1
  vAa <- 1 - s
  vaa <- 1 - (2*s)
  
  for (r in 1:nrepl) {
    
    for (i in 2:t) {
      
      # selection
      fpA <- fA[r,i-1] * (2*vAA*fA[r,i-1] + (vAa*(1-fA[r,i-1]))) / (vAA*fA[r,i-1]^2 + 2*vAa*fA[r,i-1]*(1-fA[r,i-1]) + vaa*(1-fA[r,i-1])^2)
      
      if (fpA <= 0) { fA[r,i:t] <- 0; break} # lost
      if (fpA >= 1) { fA[r,i:t] <- 1; break} # fixed
      
      # drift
      fA[r,i] <- sum(sample(x=c(0,1), size=(2*N), replace=T, prob=c((1-fpA),fpA))) / (2*N)
      
    }
    
  }
  
  u <- 0
  if ((2*N*s) > -1) u <- 1/(2*N)
  if ((2*N*s) > 1) u <- 2*s
  
  cat("Lost = ", length(which(fA[,t]==0)), "\n")
  cat("Fixed = ", length(which(fA[,t]==1)), "\t (expected = ", (u*nrepl), ")\n")
  
  return(invisible(fA));
  
}

##### :)
plotTrajectory <- function(fA, ylim=c(0,1), tlim=c(1,NA)) {
  cols <- colors()
  if (is.na(tlim[2])) tlim <- c(1,ncol(fA))
  plot(fA[1,],ylim=ylim,ty="l",xlim=tlim,col=cols[2],xlab="generations",ylab="frequency",lwd=2)
  for (i in 2:nrow(fA)) lines(fA[i,],type="l",col=cols[i+1],lwd=2)
}

plotTrajectory(simulateTrajectory(s=0.001, N=100, t=100, nrepl=100))



