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

    jc <- coda.samples(jm, jp, n.iter=12000)

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
##                    n=n.unmarked,   ## BUG FIX: 2021-08-10
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

    jc2 <- coda.samples(jm2, jp2, n.iter=12000)

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


## Doing batch 1 2021-09-13 11:35:59 
## Doing batch 2 2021-09-14 10:51:18 
## Doing batch 3 2021-09-15 10:24:43 
## Doing batch 4 2021-09-16 09:27:51 
## Doing batch 5 2021-09-17 08:05:15 
## Doing batch 6 2021-09-18 11:00:07 
## Doing batch 7 2021-09-19 09:42:16 
## Doing batch 8 2021-09-20 06:06:00 
## Doing batch 9 2021-09-21 05:58:25 
## Doing batch 10 2021-09-22 02:51:05 
## Done 2021-09-23 01:46:30 



## Doing batch 1 2021-06-09 16:06:33
## cat("Done", format(Sys.time()), "\n")
## Doing batch 2 2021-06-10 13:03:53
## Doing batch 3 2021-06-11 09:14:23
## Doing batch 4 2021-06-12 03:53:16
## Doing batch 5 2021-06-13 01:06:20
## Doing batch 6 2021-06-14 12:08:49
## Doing batch 7 2021-06-16 03:09:37
## Doing batch 8 2021-06-17 02:45:22
## Doing batch 9 2021-06-17 21:29:46
## Doing batch 10 2021-06-18 19:44:27
## > Done 2021-06-19 15:10:14


## > source("fit_case1_two-stage.R")
## Doing batch 1 2021-07-06 11:52:33 
## Doing batch 2 2021-07-07 12:21:00 
## Doing batch 3 2021-07-08 12:39:58 
## Doing batch 4 2021-07-09 13:09:19 
## Doing batch 5 2021-07-10 13:33:46 
## Doing batch 6 2021-07-11 14:05:51 
## Doing batch 7 2021-07-12 17:54:47 
## Doing batch 8 2021-07-13 19:55:52 
## Doing batch 9 2021-07-14 20:42:30 
## Doing batch 10 2021-07-16 13:29:31 
## Done 2021-07-18 04:21:47 


## Doing batch 1 2021-08-11 08:24:55 
## save(samples.stage2, file="samples_stage2.gzip")
## cat("Done", format(Sys.time()), "\n")
## Doing batch 2 2021-08-11 23:17:00 
## Doing batch 3 2021-08-12 14:03:58 
## Doing batch 4 2021-08-13 04:57:20 
## Doing batch 5 2021-08-13 19:39:11 
## Doing batch 6 2021-08-14 10:10:28 
## Doing batch 7 2021-08-15 00:50:28 
## Doing batch 8 2021-08-15 15:45:19 
## Doing batch 9 2021-08-16 06:31:28 
## Doing batch 10 2021-08-16 21:19:42 
## > > Done 2021-08-17 12:08:01 
