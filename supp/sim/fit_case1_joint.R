library(parallel)


load("sims-case1.gzip")



n.cores <- 10
(n.datasets <- length(sims.case1))
(n.batches <- ceiling(n.datasets/n.cores))
(batch.size <- n.datasets/n.batches)


cl <- makeCluster(n.cores)

## This function will analyze each dataset when called from a unique core
do.jags <- function(dataset) {

    M <- 175 ## You can change M to anything less than the M used to simulate data
    jd <- with(dataset,
               list(M=M, 
                    K = dim(u)[2],
                    J=nrow(n.all), T=ncol(n.all),
                    xlim = xlim, ylim = ylim, 
                    Area = Area, x = x, n.marked=n.marked,
                    n=n.unmarked,
                    y=y.marked[1:M,,],
                    y.cap=y.cap[1:M,],
                    u=u[1:M,,,]))

    ji <- function() {
        ## Initialize activity centers for marked guys
        si <- dataset$latent$s
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
        list(s=si[1:jd$M,,], z=matrix(1,jd$M,jd$T), 
             beta0.ED = runif(1, 0.1, 0.5), beta1.ED = 0, 
             sigmaMean=runif(1, 600, 700),
             lam0Mean=runif(1, 0.02, 0.03),
             p.cap=runif(1, 0, 0.5), alpha=0,
             eps.sd=0.1, epsilon=rep(0, jd$T) )
    }
    
    jp <- c("sigmaMean", "lam0Mean", "beta0.ED", "beta1.ED", "alpha", "N", "ED", "p.cap")

    library(rjags)

    jm <- jags.model(file="gsmr.jag", data=jd, inits=ji, n.adapt=100)

    jc <- coda.samples(jm, jp, n.iter=12000)

    return(jc)

}

samples.joint <- vector(mode="list", length=n.datasets)

for(i in 1:n.batches) {
    cat("Doing batch", i, format(Sys.time()), "\n")
    batch <- seq(1+batch.size*(i-1), length.out=batch.size)
    samples.joint[batch] <- parSapply(cl=cl, X=sims.case1[batch],
                                      FUN=do.jags)
    if(i==1)
        save(samples.joint, file="samples_joint.gzip")
    gc()
}

save(samples.joint, file="samples_joint.gzip")

cat("Done", format(Sys.time()), "\n")


## Doing batch 1 2021-09-13 11:34:17 
## Doing batch 2 2021-09-14 11:22:43 
## Doing batch 3 2021-09-15 10:39:33 
## Doing batch 4 2021-09-16 10:17:28 
## Doing batch 5 2021-09-17 12:01:25 
## Doing batch 6 2021-09-19 14:51:36 
## Doing batch 7 2021-09-22 00:47:44 
## Doing batch 8 2021-09-24 20:49:15 
## Doing batch 9 2021-09-27 16:45:55 
## Doing batch 10 2021-09-30 13:12:25 
## Done 2021-10-03 12:41:23 


## Doing batch 1 2021-06-09 16:04:55 
## Doing batch 2 2021-06-10 15:27:33 
## Doing batch 3 2021-06-11 13:14:52 
## Doing batch 4 2021-06-12 10:34:51 
## Doing batch 5 2021-06-13 09:38:19 
## Doing batch 6 2021-06-15 12:00:06 
## Doing batch 7 2021-06-17 20:23:14 
## Doing batch 8 2021-06-20 20:17:56 
## Doing batch 9 2021-06-23 17:41:55 
## Doing batch 10 2021-06-26 15:19:31 
## Done 2021-06-29 12:36:20 


## Doing batch 1 2021-07-06 11:48:07 
## Doing batch 2 2021-07-07 04:15:16 
## Doing batch 3 2021-07-07 20:54:37 
## Doing batch 4 2021-07-08 13:27:56 
## Doing batch 5 2021-07-09 08:20:28 
## Doing batch 6 2021-07-11 03:58:38 
## Doing batch 7 2021-07-13 08:42:23 
## Doing batch 8 2021-07-16 01:46:52 
## Doing batch 9 2021-07-18 23:03:08 
## Doing batch 10 2021-07-21 20:48:11 
## Done 2021-07-24 19:11:04 
