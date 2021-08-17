
sim.gsmr <- function(trapLocs,
                     buffer,
                     nPrimary,
                     nSecondary,
                     beta0, beta1,
                     alpha, epsSD,
                     lam0Mean, sigmaMean,
                     p.cap,
                     M) {
    
    nTraps <- nrow(trapLocs)
    xlim <- range(trapLocs[,1])+c(-buffer,buffer)
    ylim <- range(trapLocs[,2])+c(-buffer,buffer)
    Area <- diff(xlim)*diff(xlim) / 1e6
    psi <- N <- EN <- ED <- rep(NA, nPrimary)
    epsilon <- rep(NA, nPrimary)
    epsilon[1] <- rnorm(1, mean=0, epsSD)
    ED[1] <- max(beta0 + beta1*1 + epsilon[1], 0)
    EN[1] <- ED[1]*Area
    psi[1] <- EN[1]/M
    ## N[1] <- rpois(1, EN[1])
    lam0 <- rep(lam0Mean, nPrimary)
    sigma <- rep(sigmaMean, nPrimary)
    for(t in 2:nPrimary) {
        epsilon[t] <- rnorm(1, alpha*epsilon[t-1], epsSD)
        ED[t] <- max(beta0 + beta1*t + epsilon[t], 0)
        EN[t] <- ED[t]*Area
        psi[t] <- EN[t]/M
        ## N[t] <- rpois(1, EN[t])
    }
    s <- array(NA, c(M, 2, nPrimary))
    z <- matrix(NA, M, nPrimary)
    lambda <- dist <- array(NA, c(M, nTraps, nPrimary))
    y.counts <- array(0L, c(M,nTraps,nPrimary))
    y.marked.aug <- array(0L, c(M, nTraps, nPrimary))
    n.all <- n.unmarked <- p <- Lambda <- matrix(NA, nTraps, nPrimary)
    u.aug <- u <- array(NA, c(M, nSecondary, nPrimary, 2))
    y.cap <- matrix(NA, M, nPrimary)
    y.cap.aug <- matrix(0L, M, nPrimary)
    y.bin <- array(NA_integer_, c(M, nTraps, nSecondary, nPrimary))
    n.marked <- rep(NA, nPrimary)
    for(t in 1:nPrimary) {
        ## Activity centers
        s[,1,t] <- runif(M, xlim[1], xlim[2])
        s[,2,t] <- runif(M, ylim[1], ylim[2])
        ## Data augmentation bit
        z[,t] <- rbinom(M, 1, psi[t])
        ## Capture process
        y.cap[,t] <- rbinom(M, 1, p.cap*z[,t])
        unmarked <- ((1-y.cap[,t])*z[,t])==1
        ## Telemetry locations (one per day)
        marked <- y.cap[,t]==1
        n.marked[t] <- sum(marked)
        for(k in 1:nSecondary) {
            u[,k,t,1] <- rnorm(M, s[,1,t], sigma[t])
            u[,k,t,2] <- rnorm(M, s[,2,t], sigma[t])
        }
        for(j in 1:nTraps) {
            dist[,j,t] <- sqrt((s[,1,t]-trapLocs[j,1])^2 +
                               (s[,2,t]-trapLocs[j,2])^2)
            lambda[,j,t] <- lam0[t]*exp(-dist[,j,t]^2 /
                                        (2*sigma[t]^2))
            for(k in 1:nSecondary) {
                y.bin[,j,k,t] <- rbinom(M, 1, 1-exp(-lambda[,j,t]*z[,t]))
            }
        }
        y.marked.aug[1:n.marked[t],,t] <- apply(y.bin[marked,,,t,drop=FALSE], c(1,2), sum)
        nk <- apply(y.bin[unmarked,,,t,drop=FALSE], c(2,3), sum)
        nkb <- ifelse(nk>0, 1, 0)
        n.unmarked[,t] <- rowSums(nkb)
        y.all <- ifelse(apply(y.bin[,,,t], c(2,3), sum)>0, 1, 0)
        n.all[,t] <- rowSums(y.all)
        y.cap.aug[1:n.marked[t],t] <- y.cap[marked,t]
        u.aug[1:n.marked[t],,t,] <- u[marked,,t,]
    }

    out <- list(pars=c(beta0, beta1, alpha, epsSD,
                       lam0Mean, sigmaMean),
                latent=list(s=s,z=z,y.bin=y.bin,y.cap=y.cap),
                lam0=lam0, sigma=sigma,
                x=trapLocs, xlim=xlim, ylim=ylim,
                y.marked=y.marked.aug, y.cap=y.cap.aug,
                u=u.aug, n.unmarked=n.unmarked, n.all=n.all,
                n.marked=n.marked, Area=Area)

    return(out)

}





traps.in <- read.csv("../../data/trap_locations.csv")

traps.al <- as.matrix(traps.in[1:60,2:3])


ED <- eps <- rep(NA, 50)
eps[1] <- rnorm(1, 0, 0.1)
ED[1] <- 1.5 + 0.005 + eps[1]
for(t in 2:50) {
    eps[t] <- rnorm(1, 0.5*eps[t-1], 0.1)
    ED[t] <- 1.5 - 0.005*t + eps[t]
}

plot(ED, type='o')



set.seed(43545)

tmp <- sim.gsmr(trapLocs=traps.al, buffer=1000,
                nPrimary=20, nSecondary=14,
                beta0=1.5, beta1=-0.005,
                alpha=0.5, epsSD=0.1,
                lam0Mean=0.07, sigmaMean=400, M=500,
                p.cap=0.1)

str(tmp)

colSums(tmp$y.cap)     ## Deer with telem
colSums(tmp$latent$z)  ## Abundance over time
table(tmp$y.marked)
table(tmp$n.unmarked)


plot(function(x) (1 - 0.02*x)*tmp$Area, from=1, to=20, ylim=c(0, 100))




## Simulate the datasets

set.seed(4300)

sims.case1 <- replicate(
    n=100,
    sim.gsmr(trapLocs=traps.al, buffer=1000,
             nPrimary=20, nSecondary=14,
             beta0=1, beta1=-0.02,
             alpha=0.5, epsSD=0.1,
             lam0Mean=0.07, sigmaMean=400, M=500,
             p.cap=0.2),
    simplify=FALSE)
                        
str(sims.case1,0)

sims.case1[[1]]$n.marked
colSums(sims.case1[[1]]$latent$z)


save(sims.case1, file="sims-case1.gzip")



## Test run in JAGS

library(rjags)

tmp <- sims.case1[[1]]

M <- 175 ## You can change M to anything less than the M used to simulate data
jd <- with(tmp,
           list(M=M, 
                K = dim(u)[2],
                J=nrow(n.all), T=ncol(n.all),
                xlim = xlim, ylim = ylim, 
                Area = Area, x = x, n.marked=n.marked,
                n=n.unmarked,
                y=y.marked[1:M,,],
                y.cap=y.cap[1:M,],
                u=u[1:M,,,]))

str(jd)



ji <- function() {
    ## Initialize activity centers by putting marked guys first
    si <- tmp$latent$s
    marked.old <- tmp$latent$y.cap==1
    n.marked.old <- colSums(marked.old)
    for(t in 1:jd$T) {
        si[1:n.marked.old[t],,t] <- si[marked.old[,t],,t]
    }
    list(s=si[1:jd$M,,], z=matrix(1,jd$M,jd$T), 
         beta0.ED = rnorm(1,2,0.5), beta1.ED = 0, 
         sigmaMean=runif(1, 700, 900),
         lam0Mean=runif(1, 0.03, 0.05),
         p.cap=runif(1), alpha=0,
         eps.sd=0.1, epsilon=rep(0, jd$T) )
}

str(ji())


jp <- c("sigmaMean", "lam0Mean", "beta0.ED", "beta1.ED", "alpha", "N", "ED")


system.time({
    jm <- jags.model(file="gsmr.jag", data=jd,
                     inits=ji, n.adapt=100)
})

system.time({
    jc <- coda.samples(jm, jp, n.iter=100)
}) ## 366 it/hr with M=175 :(







