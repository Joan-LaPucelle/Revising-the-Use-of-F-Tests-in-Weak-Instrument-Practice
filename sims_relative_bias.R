rm(list=ls())

# Simulations - Relative bias of 2SLS estimates
# Ashish Patel - notes for James Lane (29/6/25)

set.seed(100)
runs = 1000 # number of simulation runs
n=50000 # sample size
J = 4 # number of instruments
beta = 1 # true causal effect
sig.v <- 1; sig.u <- 1; rho <- 0.5 # error variances and covariance
sim.sample <- 100000 # simulation draws to find the relative bias level within each simulation run 

# finding the value of l such that the relative bias is 0.1
B <- function(l){
  z_v <- mvtnorm::rmvnorm(sim.sample,mean=rep(0,J),sigma=diag(J)) # sim.sample number of draws of standard normal rvs
  B_n <- function(s){as.numeric(((t(z_v[s,])%*%z_v[s,])+(sqrt(l)*sum(z_v[s,])))/((t(z_v[s,])%*%z_v[s,])+(2*sqrt(l)*sum(z_v[s,]))+(J*l)))} # relative bias
  B_n <- sapply(1:sim.sample, B_n) # relative bias evaluated at each simulation run
  return(abs(mean(B_n))-0.1) # measures how far away the average relative bias is away from 0.1
}
system.time(l_opt <- uniroot(B, lower=0.01, upper=max(20,J), tol=1e-4))
l <- l_opt$root # this is value of l that corresponds to a relative bias of 0.1

# F test critical value based on Staiger and Stock asymptotics
lb <- qchisq(0.95,df=J,ncp=(J*l))/J # given the distribution of the F-statistic (p.214 of SW16), this is the critical value for F-test for a 5% test

# for simulation study, we choose the value of IV strength c needed for relative bias to be 0.1
c = sqrt(l*(sig.v^2))
pi = rep(c,J)/sqrt(n) # instrument effects on exposure under weak IV asymptotics

# simulation function
sim.fun <- function(s){
e = mvtnorm::rmvnorm(n,mean=c(0,0),sigma=matrix(c(sig.u^2,rho,rho,sig.v^2),2,2)) # correlated normally distributed errors
u = e[,1]; v = e[,2]; rm(e) # errors for exposure and outcome models
Z <- matrix(NA,nrow=n,ncol=J)
for (j in 1:J){Z[,j] <- rnorm(n,0,1)} # generate J instruments
X = Z%*%pi + v; X = as.vector(X) # generate exposure
Y = X*beta + u; Y = as.vector(Y) # generate outcome

# F statistic
sig_v <- as.numeric(var(X)-(t(cov(Z,X))%*%solve(var(Z))%*%cov(Z,X))) # estimate of first-stage residual variance (denoted small-sigma_v^2 in the paper)
W_pi <- n*as.numeric(t(cov(Z,X))%*%solve(var(Z))%*%cov(Z,X))/sig_v 
F <- W_pi/J # F-statistic for testing instrument relevance

# OLS estimator
beta_ols = cov(X,Y)/var(X)

# two-stage least squares (TSLS) estimator
beta_tsls = as.numeric(t(cov(Z,X))%*%solve(var(Z))%*%cov(Z,Y))/as.numeric(t(cov(Z,X))%*%solve(var(Z))%*%cov(Z,X))
sig_u <- mean((Y-X*beta_tsls)^2) # estimate of second-stage residual variance (denoted small-sigma_u^2 in the paper) 
W <- function(beta0){n*((beta_tsls-beta0)^2)*as.numeric(t(cov(Z,X))%*%solve(var(Z))%*%cov(Z,X))/sig_u} # Wald statistic under the null hypothesis H0: beta=beta0

if (s%%100==0){print(paste("percentage of simulations completed:", (s/runs)*100))}

return(list("ols"=beta_ols,"tsls"=beta_tsls,"W"=W(beta), "F"=F))
}

# run simulations
sims <- sapply(1:runs,sim.fun)

# collect simulation results
ols <- vector(,length=runs)
tsls <- vector(,length=runs)
wald <- vector(,length=runs)
fstat <- vector(,length=runs)
for (s in 1:runs){
  ols[s] <- sims[,s]$ols
  tsls[s] <- sims[,s]$tsls
  wald[s] <- sims[,s]$W
  fstat[s] <- sims[,s]$F
}

# ols and 2sls estimates
round(c(mean(ols),mean(tsls)),4)

# standard deviation of ols and 2sls estimates
round(c(sd(ols),sd(tsls)),4)

# type I error rate of conventional Wald test
mean(wald>qchisq(0.95,1))

# relative bias (this should be close to 0.1)
round(abs(mean(tsls-beta)/mean(ols-beta)),3)

# F statistic
round(mean(fstat),3)

# standard deviation of F statistic estimates
round(sd(fstat),3)

# F-test type I error rates using critical values from usual (strong instrument) asymptotics
mean(J*fstat>qchisq(0.95,J)) # note that we are pretty much always rejecting the null hypothesis that instruments are not relevant
# ... but of course, it's more informative to learn whether instruments are strong enough that the relative bias is at most 0.1. Lets test that now. 

# F-test type I error rates using critical values based on Staiger and Stock's weak instrument asymptotics
mean(fstat>lb) # this should be close to 0.05 since we've chosen c (instrument strength) such that the relative bias is exactly 0.1

fname = paste('relativebias_',n,'.RData', sep="")
save.image(fname)