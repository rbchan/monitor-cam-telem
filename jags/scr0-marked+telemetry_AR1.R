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
load("deer.RData")

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

load.module("dic") # Monitor deviance

n.adapt = 100 # Set adaptive phase
n.iter = 15000 # Set number of iterations

# Configure the model 
jm.bern.auto.does <- jags.model(file="scr0-marked+telemetry-AR1.jag",
                                data=jd.bern.auto.does,
                                inits=ji.bern.auto.does,
                                n.adapt=n.adapt)


# MCMC sampler
jc.bern.auto.does <- coda.samples(jm.bern.auto.does, jp.bern.auto.does,
                                  n.iter=n.iter)

# Exclude burnin
burnin <- 5001 + n.adapt
jc.bern.auto.does.nob <- window(jc.bern.auto.does, start=burnin)
