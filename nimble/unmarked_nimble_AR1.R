#### R Code Description - Margenau et al. () - Ecological Applications ####
### Authors: Lydia L. S. Margenau and R. B. Chandler
### Last updated: 24 February 2020
### Description: Autoregressive unmarked model for camera data for Bear Island.
###              Implementation of SCR model using Nimble MCMC sampler  

#### Load in required packages and data ####
# Packages needed
install.packages(c("nimble", "coda"), dependencies = TRUE)
library(nimble)
library(coda)

#### Load in deer unmarked data + marked prior distribution means 
load("unmarked_data.RData")

### Detection parameter priors from marked posterior means 
prior_means  # [t,(sigma, lam0)]
prior_varcov  # [t,(sigma, lam0),(sigma, lam0) ]

range <- match(rownames(prior_means), biweek) 
T <- length(range) # num of fortnights 
M <- 500 # Augmentation
K <- 14  # number of secondary sampling occasions
x # trap locations 
J <- nrow(x) # number of camera traps 

# Detection histories - detected/not detected counts of unmarked female deer 
# matrix(trap, sampling occasions)
n 

# Split into 3D array
# array(trap, secondary sampling occasion, primary time period)
n3D 

# Number of detections per camera trap during each fortnight 
nsum <- apply(n3D, c(1,3), sum)


# 3D operation matrix - array(trap, secondary sampling, primary fortnight)
oper

# Number of camera days per camera per fortnight 
K.oper <- apply(oper, c(1,3), sum)

# Boundaries + Area
xlim = range(x[,1])+c(-1000,1000)
ylim = range(x[,2])+c(-1000,1000)
Area <- (xlim[2]-xlim[1])*(ylim[2]-ylim[1])/1e6 ## sq-km

#### Model Developement ####

densityNimble <- nimbleCode({
  
  ### Prior Distributions 
  # Autoregressive estimation 
  for(t in 1:T) { 
    # Use the posterior estimates from the marked population as informative 
    # priors for the unmarked population
    # Multivariate normal distribution
    log_sig_lam0[1:2,t] ~ dmnorm(prior_means[t,1:2], cov=prior_vcov[t,1:2,1:2])
    sigma[t] <- exp(log_sig_lam0[1,t])
    lam0[t] <- exp(log_sig_lam0[2,t])
  }
  
  beta0.ED ~ dnorm(5, sd = 2.5)  # intercept of ED
  beta1.ED ~ dnorm(0, sd = 0.1)  # slope of ED     
  eps.sd ~ dunif(0, 2) ## SD of density noise
  alpha ~ dunif(-1, 1) ## AR1 correlation
  epsilon[1] ~ dnorm(0, sd=eps.sd) # density noise
  ED[1] <- max(beta0.ED + beta1.ED*1 + epsilon[1], 0) # left truncated glm  
  EN[1] <- ED[1]*Area # convert ED to EN
  psi[1] <- EN[1]/M # proportion of population in sample 
  
  for(t in 2:T){
    epsilon[t] ~ dnorm(alpha*epsilon[t-1], sd=eps.sd) # autoregression on noise 
    ED[t] <- max(beta0.ED + beta1.ED*t + epsilon[t], 0) 
    EN[t] <- ED[t]*Area
    psi[t] <- EN[t]/M
  }
  

   ###Unmarked SCR

   for(t in 1:T){ # loop over each fortnight
     for(i in 1:M){ # loop over each individual

       # probability of being a member of the population
       z[i,t] ~ dbern(psi[t])

       # activity centers probabilities assummed uniform
       s[i,1,t] ~ dunif(xlim[1], xlim[2])
       s[i,2,t] ~ dunif(ylim[1], ylim[2])

       for(j in 1:J){ # loop over each trap location

         # calculate distance from activity center to trap location
         d[i,j,t] <- sqrt((s[i,1,t]-x[j,1])^2 + (s[i,2,t]-x[j,2])^2)
         lambda[i,j,t] <- lam0[t]*exp(-d[i,j,t]^2 / (2*sigma[t]^2))
       }
     }

     for(j in 1:J) {
       Lambda[j,t] <- inprod(lambda[1:M,j,t] , z[1:M,t])
       p[j,t] <- 1-exp(-Lambda[j,t])
       nsum[j,t] ~ dbin(p[j,t], K[j,t])
     }
     N[t] <- sum(z[1:M,t])
     D[t] <- N[t]/Area
   }
 })


#### Model Intialization ####

constants <- list(M=M, K = K.oper, J=J, T=T,
                  xlim = xlim, ylim = ylim, 
                  Area = Area, x = x)

data <- list(nsum = nsum, 
             prior_means = prior_means, prior_vcov = prior_varcov)

# Initial values 
log.sig.lam0.init <- matrix(NA, 2, T)
log.sig.lam0.init[1,] <- log(800)
log.sig.lam0.init[2,] <- log(0.1)

M.x <- runif(M, xlim[1], xlim[2])
M.y <- runif(M, ylim[1], ylim[2])

si <- array(NA, c(M, 2, T))
for(i in 1:M){
  si[i, 1, ] <- rep(M.x[i], T)
  si[i, 2, ] <- rep(M.y[i], T)
}


inits <- list(s = si, z = matrix(1,M,T), epsilon = rep(0, T),  
              beta0.ED = rnorm(1,2,0.5), beta1.ED = rnorm(1,0,0.05), 
              eps.sd = runif(1,0,2), alpha = runif(1,-1,1),
              log_sig_lam0 = log.sig.lam0.init) 


params <- c("sigma", "lam0", "psi", 
            "eps.sd", "alpha", "epsilon", 
            "beta0.ED", "beta1.ED", 
            "EN","N", "ED", "D")

#Single chain 
system.time(mod.nimble <- nimbleModel(code = densityNimble,
                          constants = constants,
                          data = data,
                          inits = inits))


#### Model Evaluation #### 
  
  # Did everything initalize correctly?
  mod.nimble$initializeInfo()
  
  # create an MCMC configuration
  mcmcSCR<-configureMCMC(mod.nimble, monitors=params)
  
  # build the MCMC object and its samplers either from the model
  SCRMCMC <- buildMCMC(mcmcSCR)
  
  # compile the MCMC object
  Cmodel <- compileNimble(mod.nimble) 
  CompSCRMCMC <- compileNimble(SCRMCMC, project = mod.nimble)
  
  niter = 25000 # total number of iterations to run
  nburnin = 5000 # size of burnin
  thin = 2 # thinning rate
  
  ### MCMC run using runMCMC(), allows for direct conversion to coda.mcmc 
  samplesList <- runMCMC(CompSCRMCMC, niter = niter, nchains = 3,
                                      thin = thin, nburnin = nburnin,
                                      samplesAsCodaMCMC = TRUE,
                                      summary = TRUE)
  


  summary(samplesList)


