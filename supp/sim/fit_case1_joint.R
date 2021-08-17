library(parallel)


load("sims-case1.gzip")



(n.cores <- min(detectCores()-1, 10))
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
}

save(samples.joint, file="samples_joint.gzip")

cat("Done", format(Sys.time()), "\n")


