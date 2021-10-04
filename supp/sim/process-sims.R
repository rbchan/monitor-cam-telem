library(coda)
library(lattice)


list.files()

load("samples_stage1.gzip")
load("samples_stage2.gzip")
load("samples_joint.gzip")


ls()


plot(samples.stage1[[1]])
plot(samples.stage2[[1]])
plot(samples.joint[[1]])


## Compare encounter rate parameters

burn <- 1:1000

lam.sig.mean.st1 <- t(sapply(samples.stage1,
                             function(x) colMeans(x[-burn,])))
lam.sig.mean.st2 <- t(sapply(samples.stage2, function(x)
    colMeans(x[-burn,c("lam0Mean", "sigmaMean")])))
lam.sig.mean.joint <- t(sapply(samples.joint, function(x)
    colMeans(x[-burn,c("lam0Mean", "sigmaMean")])))

## Stage 1 estimates
hist(lam.sig.mean.st1[,1]); abline(v=0.07, col=4, lwd=3)
hist(lam.sig.mean.st1[,2]); abline(v=400, col=4, lwd=3)

## Stage 2 estimates (updated)
hist(lam.sig.mean.st2[,1]); abline(v=0.07, col=4, lwd=3)
hist(lam.sig.mean.st2[,2]); abline(v=400, col=4, lwd=3)

## Joint estimates 
hist(lam.sig.mean.joint[,1]); abline(v=0.07, col=4, lwd=3)
hist(lam.sig.mean.joint[,2]); abline(v=400, col=4, lwd=3)

## Same result seen here
plot(lam.sig.mean.st1[,1], lam.sig.mean.joint[,1])
plot(lam.sig.mean.st1[,2], lam.sig.mean.joint[,2])

## And here. 
pdf("sim-lam0sig.pdf", width=4, height=7)
op <- par(mfcol=c(2,1), mai=c(0.7, 0.7, 0.4, 0.1), mgp=c(2.2,1,0))
plot(lam.sig.mean.st2[,1], lam.sig.mean.joint[,1],
     xlab="Two-stage", ylab="Joint", main="lam0"); abline(0,1)
plot(lam.sig.mean.st2[,2], lam.sig.mean.joint[,2],
     xlab="Two-stage", ylab="Joint", main="sigma"); abline(0,1)
par(op)
dev.off()
system("gopen sim-lam0sig.pdf")




plot(lam.sig.mean.st1, type="n", xlim=c(0.06,0.08))
## points(lam.sig.mean.st2, pch=2, col=2)
points(0.07, 400, pch=16, cex=2, col="red")
arrows(lam.sig.mean.st1[,1], lam.sig.mean.st1[,2],
       lam.sig.mean.st2[,1], lam.sig.mean.st2[,2],
       col=rgb(0,0,0,0.2), length=0.08)

plot(lam.sig.mean.st2, xlim=c(0.05, 0.09), ylim=c(385, 415))
points(lam.sig.mean.joint, pch=2, col=2)




## Compare beta parameters

beta.mean.st2 <- t(sapply(samples.stage2, function(x)
    colMeans(x[-burn,c("beta0.ED", "beta1.ED")])))
beta.mean.joint <- t(sapply(samples.joint, function(x)
    colMeans(x[-burn,c("beta0.ED", "beta1.ED")])))

## Stage 2 estimates
hist(beta.mean.st2[,1]); abline(v=1, col=4, lwd=3)
hist(beta.mean.st2[,2]); abline(v=-0.02, col=4, lwd=3)

## Joint estimates 
hist(beta.mean.joint[,1]); abline(v=1.0, col=4, lwd=3)
hist(beta.mean.joint[,2]); abline(v=-0.02, col=4, lwd=3)

## And here. 
pdf("sim-betas.pdf", width=4, height=7)
op <- par(mfcol=c(2,1), mai=c(0.7, 0.7, 0.4, 0.1), mgp=c(2.2,1,0))
plot(beta.mean.st2[,1], beta.mean.joint[,1],
     xlab="Two-stage", ylab="Joint", main="beta0"); abline(0,1)
plot(beta.mean.st2[,2], beta.mean.joint[,2],
     xlab="Two-stage", ylab="Joint", main="beta1"); abline(0,1)
par(op)     
dev.off()
system("gopen sim-betas.pdf")

plot(beta.mean.st2, xlim=c(0.5, 1.5), ylim=c(-0.05, 0))
points(beta.mean.joint, pch=2, col=2)

## Outlier
which( beta.mean.st2[,2]< -0.15)

plot(samples.stage2[[67]][,c("beta0.ED","beta1.ED")])
plot(window(samples.joint[[67]][,c("beta0.ED","beta1.ED")],
            start=101, end=2100))


## Compare N(t)

N.mean.st2 <- t(sapply(samples.stage2, function(x)
    colMeans(x[,grep("N\\[", varnames(x))])))
N.low.st2 <- t(sapply(samples.stage2, function(x)
    apply(x[,grep("N\\[", varnames(x))], 2, quantile, prob=0.025)))
N.upp.st2 <- t(sapply(samples.stage2, function(x)
    apply(x[,grep("N\\[", varnames(x))], 2, quantile, prob=0.975)))
N.mean.joint <- t(sapply(samples.joint, function(x)
    colMeans(x[,grep("N\\[", varnames(x))])))
N.low.joint <- t(sapply(samples.joint, function(x)
    apply(x[,grep("N\\[", varnames(x))], 2, quantile, prob=0.025)))
N.upp.joint <- t(sapply(samples.joint, function(x)
    apply(x[,grep("N\\[", varnames(x))], 2, quantile, prob=0.975)))

N.mean.st2.df <- data.frame(N=as.vector(N.mean.st2),
                            low=as.vector(N.low.st2),
                            upp=as.vector(N.upp.st2),
                            year=rep(1:20, each=nrow(N.mean.st2)),
                            sim=as.character(rep(1:100, times=ncol(N.mean.st2))))
N.mean.joint.df <- data.frame(N=as.vector(N.mean.joint),
                              low=as.vector(N.low.joint),
                              upp=as.vector(N.upp.joint),
                              year=rep(1:20, each=nrow(N.mean.joint)),
                              sim=as.character(rep(1:100, times=ncol(N.mean.joint))))

## Stage 2 estimates
histogram(~N|year, N.mean.st2.df)

## Joint estimates 
histogram(~N|year, N.mean.joint.df)

plot(N.mean.st2[,1], N.mean.joint[,1])



## Compare to data-generating (actual) values of abundance

load("sims-case1.gzip")



N.true <- t(sapply(sims.case1, function(x) colSums(x$latent$z)))

N.true.df <- data.frame(N=as.vector(N.true),
                        low=NA,
                        upp=NA,
                        year=rep(1:20, each=nrow(N.true)),
                        sim=as.character(rep(1:100, times=ncol(N.true))))

N.mean <- cbind(rbind(N.mean.st2.df, N.mean.joint.df, N.true.df),
                method=rep(c("two-stage", "joint", "actual"),
                           each=nrow(N.mean.st2.df)))


pdf("sim-Nt.pdf", width=7, height=8)
xyplot(N ~ year | sim, N.mean, group=method, layout=c(5,20),
       xlab="Time",
       strip=FALSE, auto.key=list(draw=TRUE, points=FALSE, lines=TRUE),
       type="l", panel=function(...) {
           panel.xyplot(...)
##           panel.polygon()
       },
       as.table=TRUE)
dev.off()
system("gopen sim-Nt.pdf")



## Summarize sim data


N <- sapply(sims.case1, function(x) colSums(x$latent$z))
n.marked <- sapply(sims.case1, function(x) x$n.marked)

pdf("sim-N-n.pdf", width=7, height=7)
par(mai=c(0.9, 0.9, 0.1, 0.1))
matplot(N, type="l", col="darkblue", ylim=c(0, 120), lty=1,
        xlab="Time", ylab="Individuals", cex.lab=1.5)
matlines(n.marked, type="l", col="lightblue", lty=1)
legend(14, 120, c("All", "Marked"), lty=1, lwd=2,
       col=c("darkblue", "lightblue"), cex=1.3)
dev.off()
system("gopen sim-N-n.pdf")

           
mean(n.marked/N)

range(n.marked)
