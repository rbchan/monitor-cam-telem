library(parallel)


load("sims-case1.gzip")



n.cores <- 10
(n.datasets <- length(sims.case1))
(n.batches <- ceiling(n.datasets/n.cores))
(batch.size <- n.datasets/n.batches)


cl <- makeCluster(n.cores)

## This function will fit stage 1 to each dataset when called from a unique core
do.jags.stage1 <- function(dataset) {

    jd <- with(dataset,
               list(K = dim(u)[2],
                    J=nrow(n.all), T=ncol(n.all),
                    xlim = xlim, ylim = ylim, 
                    x = x, n.marked=n.marked,
                    y=y.marked[1:max(dataset$n.marked),,],
                    u=u[1:max(dataset$n.marked),,,]))

    ji <- function() {
        ## Initialize activity centers for marked guys
        ## si <- dataset$latent$s
        si <- array(NA, c(max(jd$n.marked), 2, jd$T))
        for(t in 1:jd$T) {
            for(i in 1:dataset$n.marked[t]) {
                trps <- jd$y[i,,t]>0
                if(any(trps))
                    si[i,,t] <- colMeans(jd$x[trps,,drop=FALSE])
                else {
                    si[i,,t] <- colMeans(jd$u[i,,t,])
                    ## Sometimes u is outside S
                    si[i,1,t] <- max(si[i,1,t], jd$xlim[1])
                    si[i,1,t] <- min(si[i,1,t], jd$xlim[2])
                    si[i,2,t] <- max(si[i,2,t], jd$ylim[1])
                    si[i,2,t] <- min(si[i,2,t], jd$ylim[2])
                }
            }
        }
        list(s=si[1:max(jd$n.marked),,], 
             sigmaMean=runif(1, 600, 700),
             lam0Mean=runif(1, 0.02, 0.03))
    }
    
    jp <- c("sigmaMean", "lam0Mean")

    library(rjags)

    jm <- jags.model(file="gsmr-stage1.jag", data=jd, inits=ji, n.adapt=100)

    jc <- coda.samples(jm, jp, n.iter=2000) ##12000)

    return(jc)

}



## This function will fit stage 2 to each dataset when called from a unique core
do.jags.stage2 <- function(dataset) {

    M <- 175 ## You can change M to anything less than the M used to simulate data
    jd2 <- with(dataset,
               list(M=M, 
                    K = dim(u)[2],
                    J=nrow(n.all), T=ncol(n.all),
                    xlim = xlim, ylim = ylim, 
                    Area = Area, x = x, 
                    n=n.all,
                    mean.log.sigma.lam0=mean.log.sigma.lam0,
                    vcov.log.sigma.lam0=vcov.log.sigma.lam0))

    ji2 <- function() {
        ## Initialize activity centers for marked guys
        si <- dataset$latent$s
        list(s=si[1:jd2$M,,], z=matrix(1,jd2$M,jd2$T), 
             beta0.ED = runif(1, 0.1, 0.5), beta1.ED = 0, 
             log.sigma.lam0=c(log(500), log(0.05)),
             alpha=0,
             eps.sd=0.1, epsilon=rep(0, jd2$T) )
    }
    
    jp2 <- c("sigmaMean", "lam0Mean", "beta0.ED", "beta1.ED", "alpha", "N", "ED")

    library(rjags)

    jm2 <- jags.model(file="gsmr-stage2.jag", data=jd2, inits=ji2, n.adapt=100)

    jc2 <- coda.samples(jm2, jp2, n.iter=2000) ##12000)

    return(jc2)

}




## Lists to hold posterior samples for each dataset
samples.stage1 <- vector(mode="list", length=n.datasets) 
samples.stage2 <- vector(mode="list", length=n.datasets)



## Loop over batches
for(i in 1:n.batches) {
    library(coda)
    cat("Doing batch", i, format(Sys.time()), "\n")
    batch <- seq(1+batch.size*(i-1), length.out=batch.size)
    samples.stage1[batch] <- parSapply(cl=cl, X=sims.case1[batch],
                                       FUN=do.jags.stage1)
    if(i == 1)
        save(samples.stage1, file="samples_stage1.gzip")
    for(j in batch) {
        samps <- log(as.matrix(samples.stage1[[j]]))
        ## Next line is a bug fix. Order was wrong before 2021-07-06
        samps <- samps[,c("sigmaMean", "lam0Mean")]
        xbar <- colMeans(samps)
        xvar <- var(samps)
        sims.case1[[j]]$mean.log.sigma.lam0 <- xbar
        sims.case1[[j]]$vcov.log.sigma.lam0 <- xvar
    }
    samples.stage2[batch] <- parSapply(cl=cl, X=sims.case1[batch],
                                       FUN=do.jags.stage2)
    if(i==1)
        save(samples.stage2, file="samples_stage2.gzip")
    gc()
}

save(samples.stage1, file="samples_stage1.gzip")
save(samples.stage2, file="samples_stage2.gzip")

cat("Done", format(Sys.time()), "\n")



