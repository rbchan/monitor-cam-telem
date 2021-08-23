#### R Code Description - Margenau et al. () - Ecological Applications ####
### Authors: Lydia L. S. Margenau and R. B. Chandler
### Last updated: 13 Nov 2020
### Description: Model development for biweekly estimates for a marked subset  
###              of the female deer population. An autoregressive model is 
###              implemented on the detection parameters. 
###              Model and data are formatted for MCMC sampling in NIMBLE

#### Load in required packages and data ####
# Packages needed
install.packages(c("nimble", "coda"), dependencies = TRUE)
library(nimble)
library(coda)

# Data for Bear Island camera trapping array 
load("marked_data.RData")

#### Data Description #### 

## Detection histories (present/absent) for each marked individual
## array[ID, trap, secondary occasion, primary period]
histories4D

# number of marked individuals
N <- dim(histories4D)[1]

# Trap locations for the specific camera grid   
# data.frame[UTME, UTMN]
x 

# Set the limits of study area to the camera grid with 1 km buffer  
xlim = range(x[,1])+c(-1000,1000)
ylim = range(x[,2])+c(-1000,1000)

# number of secondary sampling occasions, fortnight periods 
K <- 14

# 4-dimensional telemetry array for marked deer
# array[Deer ID x telemetry occasion x coordinates x primary period]
telemetry.deer

# Number of telemetry locations for each marked individual during each 
# primary period 
# matrix(ID, primary period)
nTelemLocs 
rowSums(nTelemLocs) # total telemetry locations per individual

# Number of primary occasions in which telemetry data is available for each 
# individual 
t.length

# The starting and ending primary period for the model based on starting 
# primary periods with marked deer 
start # Delay in marking deer on the camera grid until 9th fortnight 
end 
range <- start:end 

# Conversion of a ragged list to a matrix
# Each row contains the sequential fortnight periods an individual deer 
# was oavailable on the camera grid for sampling 
# The first column indicates when individual entered the camera array 
t.in

# Camera operational status matrix
# 1 indicates cameras was operational, 0 indicates the camera wasn't available
# array[trap, secondary occasion, primary period]
oper3D 



#### Model development ####

markedNimble <- nimbleCode({
  
  ### Prior Distributions for autoregressive terms 
  # Detection parameter- lam0
  beta0.loglam0 ~ dnorm(log(0.1), sd = 0.5) # mean (over time) lam0 on log scale  
  beta1.loglam0 ~ dnorm(0, sd = 0.1) 
  alpha.lam0  ~ dunif(-1, 1) ## AR1 correlation for lam0
  eps.sd.lam0 ~ dunif(0, 2) ## SD of lam0 noise
  epsilon.lam0[1] ~ dnorm(0, sd = eps.sd.lam0) # lam0 noise
  
    
  # spatial scale parameter - sigma 
  beta0.logsigma ~ dnorm(log(350), sd = 0.5) # mean (over time) lam0 on log scale  
  beta1.logsigma ~ dnorm(0, sd = 0.1)
  eps.sd.sigma ~ dunif(0, 2) ## SD of sigma noise
  alpha.sigma ~ dunif(-1, 1) ## AR1 correlation for sigma 
  epsilon.sigma[1] ~ dnorm(0, sd = eps.sd.sigma) # sigma noise
  
  
  # Starting values 
  loglam0[1] <- beta0.loglam0 + beta1.loglam0 * 1 + epsilon.lam0[1]
  lam0[1] <- exp(loglam0[1])
  logsigma[1] <- beta0.logsigma + beta1.logsigma * 1 + epsilon.sigma[1]
  sigma[1] <- exp(logsigma[1])
  
  #Autoregression 
  for(t in 2:length(T)){
    # autoregression on noise for lam0
    epsilon.lam0[t] ~ dnorm(alpha.lam0 * epsilon.lam0[t-1], sd = eps.sd.lam0) 
    loglam0[t] <- beta0.loglam0 + beta1.loglam0 * t + epsilon.lam0[t]
    lam0[t] <- exp(loglam0[t])
    # autoregression on noise for sigma
    epsilon.sigma[t] ~ dnorm(alpha.sigma * epsilon.sigma[t-1], sd = eps.sd.sigma)  
    logsigma[t] <- beta0.logsigma + beta1.logsigma * t + epsilon.sigma[t]
    sigma[t] <- exp(logsigma[t])
  }
  
  
  # SCR model 
  for(i in 1:n0) { # loop over ear-tagged deer
    #for(t in first[i]:last[i]) { # loop over each time interval 
    for(t in t.in[i,1]:t.in[i,t.length[i]]) { # loop over every occurence of deer on camera grid
      s[i,1,t] ~ dunif(xlim[1], xlim[2])  # activity centers probabilities assummed uniform  
      s[i,2,t] ~ dunif(ylim[1], ylim[2])
      # Model telemetry data conditional on activity centers 
      for(k in 1:nTelemLocs[i,t]) { # for each telemetry location 
        u[i,k,1,t] ~ dnorm(s[i,1,t], 1/(sigma[(start-1)]^2)) # x coord inform activity centers 
        u[i,k,2,t] ~ dnorm(s[i,2,t], 1/(sigma[(start-1)]^2)) # y coord 
      }
      for(j in 1:nTraps) { # loop over each trap location 
        # calculate distance from activity center to trap location  
        d[i,j,t] <- sqrt((s[i,1,t]-x[j,1])^2 + (s[i,2,t]-x[j,2])^2)
        for(k in 1:nOccasions) { # loop over each sampling occassion 
          # calculate encounter rate as a distance decay fx with operational matrix 
          lambda[i,j,k,t] <- lam0[(start-1)]*exp(-d[i,j,t]^2 / (2*sigma[(start-1)]^2))*oper[j,k,t] 
          # convert lambda to detection prob, Pr(at least one detection)  
          psi[i,j,k,t] <- 1-exp(-lambda[i,j,k,t]) 
          # This implies a Poisson distribution for independent counts
          # convert to binary counts using bernoulli distribution
          y[i,j,k,t] ~ dbern(psi[i,j,k,t]) 
        }
      }
    }
  }
  
})


#### Model Intialization ####

constants <- list(t.in = t.in, t.length = t.length,
                  xlim = xlim, ylim = ylim, 
                  start = start, end = end,
                  nTraps = nrow(x), T = T,
                  n0=N, nOccasions = K,  
                  x = x, 
                  nTelemLocs = nTelemLocs,
                  oper = oper3D)

data <- list(y=histories4D,  
             u=telemetry.deer 
            )

# Update activity centers 
si <- array(NA, c(N, 2,78))

for(i in 1:N) {
    si[i,1,] <- rep(runif(1, xlim[1], xlim[2]), 78) 
    si[i,2,] <- rep(runif(1, ylim[1], ylim[2]), 78)
}



# Initial values 
inits <- list(epsilon.lam0 = rep(0, length(T)), 
              epsilon.sigma = rep(0, length(T)),  
              beta0.loglam0 = runif(1,log(0.05),log(0.6)), 
              beta1.loglam0 = rnorm(1,0,0.2), 
              eps.sd.lam0 = runif(1,0,2), 
              alpha.lam0 = runif(1,-1,1),
              beta0.logsigma = runif(1,log(100),log(850)), 
              beta1.logsigma = rnorm(1,0,0.2), 
              eps.sd.sigma = runif(1,0,2), 
              alpha.sigma = runif(1,-1,1), 
              s = si
              ) 

params <- c("sigma", "lam0",  
            "beta0.loglam0", "beta1.loglam0", 
            "beta0.logsigma", "beta1.logsigma", 
            "alpha.lam0", "epsilon.lam0",
            "alpha.sigma", "epsilon.sigma", 
            "eps.sd.lam0", "eps.sd.sigma"
            )


#### Model Evaluation #### 

system.time(mod.nimble <- nimbleModel(code = markedNimble,
                          constants = constants,
                          data = data,
                          inits = inits))

# Did everything initalize correctly?
mod.nimble$initializeInfo()

# create an MCMC configuration
mcmcSCR<-configureMCMC(mod.nimble, monitors=params)
# build the MCMC object and its samplers either from the model
SCRMCMC <- buildMCMC(mcmcSCR)
# compile the MCMC object
Cmodel <- compileNimble(mod.nimble) 
CompSCRMCMC <- compileNimble(SCRMCMC, project = mod.nimble)

niter = 15000 # total number of iterations to run
nchains = 3 # number of chains to run
nburnin = 5000 # size of burnin

### MCMC run using runMCMC(), allows for direct conversion to coda.mcmc 
system.time(samplesList <- runMCMC(CompSCRMCMC, niter = niter, nchains = nchains,
                        nburnin = nburnin,
                        samplesAsCodaMCMC = TRUE,
                        summary = TRUE))


modelSummary <- samplesList$summary


