#### R Code Description - Margenau et al. () - Ecological Applications ####
### Authors: Lydia L. S. Margenau and R. B. Chandler
### Last updated: 30 April 2019
### Description: Model development for biweekly estimates for a marked subset  
###              of the female deer population. An autoregressive model is 
###              implemented on the detection parameters. 
###              Model and data are formatted for MCMC sampling in JAGS

#### Load in required packages and data ####
# Packages needed
install.packages(c("rjags", "coda"), dependencies = TRUE)
library(rjags)
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

#### Model Development ####
cat("
  # Known number of marked individuals 
  # AR1 joint model for camera and telemetry data 
  model {
    
    ### Prior Distributions for autoregressive terms 
    # Detection parameter - lam0
    beta0.loglam0 ~ dnorm(log(0.1), 1/(0.5^2)) # mean (over time) lam0 on log scale  
    beta1.loglam0 ~ dnorm(0, 0.1) 
    alpha.lam0  ~ dunif(-1, 1) ## AR1 correlation for lam0
    eps.sd.lam0 ~ dunif(0, 2) ## SD of lam0 noise
    epsilon.lam0[start] ~ dnorm(0, 1/(eps.sd.lam0^2)) # lam0 noise
    
    
    # spatial scale parameter - sigma 
    beta0.logsigma ~ dnorm(log(350), 1/(0.5^2)) # mean (over time) lam0 on log scale  
    beta1.logsigma ~ dnorm(0, 0.1)
    eps.sd.sigma ~ dunif(0, 2) ## SD of sigma noise
    alpha.sigma ~ dunif(-1, 1) ## AR1 correlation for sigma 
    epsilon.sigma[start] ~ dnorm(0, 1/(eps.sd.sigma^2)) # sigma noise
    
    
    # Starting values 
    loglam0[start] <- beta0.loglam0 + beta1.loglam0 * start + epsilon.lam0[start]
    lam0[start] <- exp(loglam0[start])
    logsigma[start] <- beta0.logsigma + beta1.logsigma * start + epsilon.sigma[start]
    sigma[start] <- exp(logsigma[start])
    
    #Autoregression 
    for(t in (start+1):end){
      # autoregression on noise for lam0
      epsilon.lam0[t] ~ dnorm(alpha.lam0 * epsilon.lam0[t-1], 1/(eps.sd.lam0^2))
      loglam0[t] <- beta0.loglam0 + beta1.loglam0 * t + epsilon.lam0[t]
      lam0[t] <- exp(loglam0[t])
      # autoregression on noise for sigma
      epsilon.sigma[t] ~ dnorm(alpha.sigma * epsilon.sigma[t-1], 1/(eps.sd.sigma^2))  
      logsigma[t] <- beta0.logsigma + beta1.logsigma * t + epsilon.sigma[t]
      sigma[t] <- exp(logsigma[t])
    }
    
    
    
    
    # SCR model 
    for(i in 1:n0) { # loop over ear-tagged deer
      for(t in t.in[i,1]:t.in[i, t.length[i]]) { # loop over every occurrence of deer on camera grid
        s[i,1,t] ~ dunif(xlim[1], xlim[2])  # activity centers probabilities, assumed uniform  
        s[i,2,t] ~ dunif(ylim[1], ylim[2])
        # Model telemetry data conditional on activity centers 
        for(k in 1:nTelemLocs[i,t]) { # for each telemetry location 
          u[i,k,1,t] ~ dnorm(s[i,1,t], 1/(sigma[t]^2)) # x coord 
          u[i,k,2,t] ~ dnorm(s[i,2,t], 1/(sigma[t]^2)) # y coord 
        }
        for(j in 1:nTraps) { # loop over each trap location 
          # calculate distance from activity center to trap location  
          d[i,j,t] <- sqrt((s[i,1,t]-x[j,1])^2 + (s[i,2,t]-x[j,2])^2)
          for(k in 1:nOccasions) { # loop over each secondary sampling occasion 
            # calculate encounter rate as a distance decay fx with camera operational matrix 
            lambda[i,j,k,t] <- lam0[t]*exp(-d[i,j,t]^2 / (2*sigma[t]^2)) * oper[j,k,t] 
            # convert lambda to detection prob, Pr(at least one detection)  
            psi[i,j,k,t] <- 1-exp(-lambda[i,j,k,t]) 
            # This implies a Poisson distribution for independent counts
            y[i,j,k,t] ~ dbern(psi[i,j,k,t]) # convert to binary counts using Bernoulli distribution
          }
        }
      }
    }
  }", file = "marked+telemetry-AR1.jag")



#### Model Initialization ####
jd.bern.auto.does <- list(y=histories4D, x=x, u=telemetry.deer,
                          xlim = xlim, ylim = ylim,
                          t.in = t.in, t.length = t.length, 
                          start = start, end = end,
                          nTraps = nrow(x), nTelemLocs = nTelemLocs,
                          n0=N, nOccasions = K, 
                          oper = oper3D)


# Initial values 
ji.bern.auto.does <- function()
  list(beta0.loglam0 = runif(1,log(0.05),log(0.6)), 
       beta1.loglam0 = rnorm(1,0,0.2), 
       eps.sd.lam0 = runif(1,0,2), 
       alpha.lam0 = runif(1,-1,1),
       
       beta0.logsigma = runif(1,log(100),log(850)), 
       beta1.logsigma = rnorm(1,0,0.2), 
       eps.sd.sigma = runif(1,0,2), 
       alpha.sigma = runif(1,-1,1)
  )



# Variables to monitor 
jp.bern.auto.does <- c("sigma", "lam0",
                       "beta0.loglam0", "beta1.loglam0", 
                       "beta0.logsigma", "beta1.logsigma", 
                       "alpha.lam0", "epsilon.lam0",
                       "alpha.sigma", "epsilon.sigma", 
                       "eps.sd.lam0", "eps.sd.sigma",
                       "deviance")


#### Model Evaluation ####

load.module("dic") # Monitor deviance

n.adapt = 100 # Set adaptive phase
n.iter = 15000 # Set number of iterations

# Configure the model 
jm.bern.auto.does <- jags.model(file="marked+telemetry-AR1.jag",
                                data=jd.bern.auto.does,
                                inits=ji.bern.auto.does,
                                n.adapt=n.adapt)


# MCMC sampler
jc.bern.auto.does <- coda.samples(jm.bern.auto.does, jp.bern.auto.does,
                                  n.iter=n.iter)

# Exclude burnin
burnin <- 5001 + n.adapt
jc.bern.auto.does.nob <- window(jc.bern.auto.does, start=burnin)
