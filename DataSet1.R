library(trust)
library(truncnorm)
library(stats4)
library(mvtnorm)
library(xtable)
library(plot3D)
library(survival)

data <- c(3, 5, 6, 7, 8, 9,10, 12,15, 18, 19, 20, 22, 25, 28, 30, 40)
cens.data <- c(10,15,45)

n = length(data)
m = length(cens.data)

########################################################################################
#                             Trust Function Optimization (MLE)                        #       
########################################################################################

objfun <- function(parameters){
  alpha <- parameters[1]
  beta <- parameters[2]
  kappa <- parameters[3]
  
  z <- 1/alpha*(cens.data^(1-kappa)/sqrt(beta) - sqrt(beta)/cens.data^kappa)
  zeta <- 1/alpha*(cens.data^(1-kappa)/sqrt(beta) + sqrt(beta)/cens.data^kappa)
  S <- 1 - pnorm(z)
  S1 <- 1/alpha*z*dnorm(z)
  S2 <- 1/(2*beta) *zeta*dnorm(z)
  S3 <- log(cens.data)*z*dnorm(z)
  S11 <- 1/alpha^2 *z*dnorm(z)*(z^2-2)
  S12 <- 1/(2*alpha*beta) *zeta*dnorm(z)*(z^2-1)
  S13 <- 1/alpha *z*log(cens.data)*dnorm(z)*(z^2-1) 
  S22 <- 1/(4*beta^2) *dnorm(z)* (z*zeta^2 - 2*zeta - z)
  S23 <- 1/(2*beta) *zeta*log(cens.data)*dnorm(z)*(z^2-1)
  S33 <- (log(cens.data))^2 *z*dnorm(z)*(z^2-1)
  
  loglikelihood <- -n*log(alpha) - n/2*log(beta) - kappa*sum(log(data)) + sum(log(1-kappa+beta*kappa/data)) - 
                    1/(2*alpha^2*beta)*sum((data-beta)^2/data^(2*kappa)) + sum(log(S))
  
  g1.comp <- -n/alpha + 1/(alpha^3*beta)*sum((data-beta)^2/data^(2*kappa))
  g1.cens <- sum(S1/S)
  
  g2.comp <- -n/(2*beta) + kappa*sum(1/(data*(1-kappa)+beta*kappa)) - 1/(2*alpha^2)*sum(data^(-2*kappa)) + 1/(2*alpha^2*beta^2)*sum(data^(2-2*kappa))
  g2.cens <- sum(S2/S)
  
  g3.comp <- -sum(log(data)) + sum((beta-data)/(data*(1-kappa)+beta*kappa)) + 1/(alpha^2*beta)*sum(data^(-2*kappa)*(data-beta)^2*log(data))
  g3.cens <- sum(S3/S)
  
  h11.comp <- n/alpha^2 - 3/(alpha^4*beta)*sum((data-beta)^2/data^(2*kappa))
  h11.cens <- sum( (S*S11 - S1^2)/S^2 )
  h11 <- h11.comp + h11.cens
  
  h22.comp <- n/(2*beta^2) - kappa^2*sum((data*(1-kappa)+beta*kappa)^(-2)) - 1/(alpha^2*beta^3)*sum(data^(2-2*kappa))
  h22.cens <- sum( (S*S22 - S2^2)/S^2 )
  h22 <- h22.comp + h22.cens 
  
  h33.comp <- -sum((beta-data)^2/(data*(1-kappa)+beta*kappa)^2) - 2/(alpha^2*beta)*sum((log(data)*(data-beta))^2/data^(2*kappa))
  h33.cens <- sum( (S*S33 - S3^2)/S^2 )
  h33 <- h33.comp + h33.cens
  
  h12.comp <- 1/alpha^3*sum(data^(-2*kappa)) - 1/(alpha^3*beta^2)*sum(data^(2-2*kappa))
  h12.cens <- sum( (S*S12 - S1*S2)/S^2 )
  h12 <- h12.comp + h12.cens 
  
  h13.comp <- -2/(alpha^3*beta)*sum(log(data)*(data-beta)^2/data^(2*kappa))
  h13.cens <- sum( (S*S13 - S1*S3)/S^2 )
  h13 <- h13.comp + h13.cens
  
  h23.comp <- sum(data/(data*(1-kappa)+beta*kappa)^2) + 1/(alpha^2*beta^2)*sum(log(data)*(beta^2-data^2)/data^(2*kappa)) 
  h23.cens <- sum( (S*S23 - S2*S3)/S^2 )
  h23 <- h23.comp + h23.cens
  
  g <- c(g1.comp+g1.cens, g2.comp+g2.cens, g3.comp+g3.cens)
  h <- rbind(c(h11,h12,h13),c(h12,h22,h23),c(h13,h23,h33))
  list(value = loglikelihood, gradient = g, hessian = h)
}
optimization <- trust(objfun,parinit=c(.9,17,.3),rinit=5,rmax=1000,blather=F,minimize=F)
(alpha.hat <- optimization$argument[1])
(beta.hat <- optimization$argument[2])
(kappa.hat <- optimization$argument[3])
(cov <- chol2inv(chol(-optimization$hessian)))

# optim.objfun <- function(alpha,beta,kappa){
#   -( -n*log(alpha) - n/2*log(beta) - kappa*sum(log(data)) + sum(log(1-kappa+beta*kappa/data)) - 
#        1/(2*alpha^2*beta)*sum((data-beta)^2/data^(2*kappa)) + 
#        sum(log(1 - pnorm(1/alpha*(cens.data^(1-kappa)/sqrt(beta) - sqrt(beta)/cens.data^kappa)))) )
# }
# optim.optimization <- mle(optim.objfun,start=list(alpha=0.8,beta=18,kappa=0.5),method="L-BFGS-B",lower=c(0.2, 5,.2),upper=c(2, 200,.8))


#hyperparameters
a0 <- 10 #must be >4 for Var(inv gamma) to exist
a1 <- 19
b0 <- 10 #must be >4 for Var(inv gamma) to exist
b1 <- 0.083
d0 <- 1
d1 <- 1

########################################################################################
#                                Gibbs Sampling                                        #
########################################################################################

#initiate parameters
nsim <- 20000
burn <- 5000 #burn-in period
alpha2 <- c(alpha.hat^2,rep(0,length=nsim)) #starting value for alpha 
beta <- c(beta.hat,rep(0,length=nsim)) #starting value for beta
kappa <- c(kappa.hat,rep(0,length=nsim)) #starting value for k
latent <- matrix(0,ncol=m,nrow=nsim) #initiate latent vars
kappa.count <- rep(1,length=nsim) #to check MH acceptance rate
beta.count <- rep(1,length=nsim) #to check MH acceptance rate

# a refers to alpha^2
log.cond.post.kappa <- function(a,b,k,latent.dat){
  - k*(sum(log(data))+sum(log(latent.dat))) + sum(log(1-k+b*k/data)) + sum(log(1-k+b*k/latent.dat)) +
    (d0-1)*log(k) + (d1-1)*log(1-k) - 1/(2*a)*( b*(a0/a1+sum(data^(-2*k))+sum(latent.dat^(-2*k))) + 1/b*(sum(data^(2-2*k))+sum(latent.dat^(2-2*k))) -
                                                  2*(sum(data^(1-2*k))+sum(latent.dat^(1-2*k))) )
}

log.cond.post.beta <- function(a,b,k,latent.dat){
  -((b0-a0+n+m)/2+1)*log(b) - b0/b1/2/b + sum(log(1-k+b*k/data)) + sum(log(1-k+b*k/latent.dat)) -
    1/(2*a)*( b*(a0/a1+sum(data^(-2*k))+sum(latent.dat^(-2*k))) + 1/b*(sum(data^(2-2*k))+sum(latent.dat^(2-2*k))) )
}

#mean for MH
kappa.meanfunc <- function(k){
  tan(pi*(k-.5))
}

# function - convert std norm into BS
# a refers to alpha^2
rootsolve.func <- function(lat.dat,x,a,b,k){
  sqrt(a)*sqrt(b)*lat.dat*x^k - x + b
}

start.time <- Sys.time()
beta.hess <- optimization$hessian[2,2] 
trans.beta.hess <- beta.hat^2*beta.hess - b0/2/b1/beta.hat - a0*beta.hat/2/a1/alpha.hat^2
trans.beta.sd <- 1/sqrt(-trans.beta.hess)

for (iter in 1:nsim){
  
  #Step 1: Sample Latent Variables from truncated normal  
  censor.thresh <- 1/sqrt(alpha2[iter])*(cens.data^(1-kappa[iter])/sqrt(beta[iter])-sqrt(beta[iter])/cens.data^kappa[iter])
  new.latent <- rtruncnorm(m,a=censor.thresh,b=Inf,mean=0,sd=1)
  latent[iter,] <- sapply(new.latent, function(new.latent) uniroot(rootsolve.func,lower=1E-3,upper=1E10,extendInt="downX",
                    lat.dat=new.latent,a=alpha2[iter],b=beta[iter],k=kappa[iter])$root)  
  
  #Step 2: Sample kappa with RW 
  
  tune <- 15
  while (kappa[iter+1]==0 | kappa[iter+1]==0.995){
    kappa[iter+1] <- rbeta(1,tune*kappa[iter], tune*(1-kappa[iter]))
  }
  log.current.posterior <- log.cond.post.kappa(alpha2[iter],beta[iter],kappa[iter],latent[iter,])
  log.proposed.posterior <- log.cond.post.kappa(alpha2[iter],beta[iter],kappa[iter+1],latent[iter,])
  log.current.proposal <- dbeta(kappa[iter],tune*kappa[iter+1], tune*(1-kappa[iter+1]),log=T) 
  log.proposed.proposal <- dbeta(kappa[iter+1],tune*kappa[iter], tune*(1-kappa[iter]),log=T)
  log.MH.ratio <- log.proposed.posterior - log.current.posterior + log.current.proposal - log.proposed.proposal
  epsilon <- min(c(1,exp(log.MH.ratio)))
  u <- runif(1)
  if(u > epsilon) {
    kappa[iter+1] <- kappa[iter]
    kappa.count[iter] <- 0
  }
  
  #Step 3: Sample beta with RW 
  
  tune.beta <- 2.4
  new.beta <- rnorm(1,mean=log(beta[iter]),sd=tune.beta*trans.beta.sd)
  beta[iter+1] <- exp(new.beta)
  log.current.posterior <- log.cond.post.beta(alpha2[iter],beta[iter],kappa[iter+1],latent[iter,])
  log.proposed.posterior <- log.cond.post.beta(alpha2[iter],beta[iter+1],kappa[iter+1],latent[iter,])
  log.MH.ratio <- log.proposed.posterior - log.current.posterior + log(beta[iter+1]) - log(beta[iter])    
  epsilon <- min(c(1,exp(log.MH.ratio)))
  u <- runif(1)
  if(u > epsilon) {
    beta[iter+1] <- beta[iter]
    beta.count[iter] <- 0
  }
  
  #Step 4: Sample alpha^2 from Conditional Inverse-Gamma Posterior  
  alpha2[iter+1] <- 1/rgamma(1, shape = (a0+n+m)/2, 
                             scale = 1/(0.5*(1/beta[iter+1]*(sum(data^(2-2*kappa[iter+1]))+sum(latent[iter,]^(2-2*kappa[iter+1]))) + 
                                               beta[iter+1]*(a0/a1+sum(data^(-2*kappa[iter+1]))+sum(latent[iter,]^(-2*kappa[iter+1]))) -
                                               2*(sum(data^(1-2*kappa[iter+1]))+sum(latent[iter,]^(1-2*kappa[iter+1]))))))
}
end.time <- Sys.time()
end.time - start.time

#Results
alpha <- sqrt(alpha2)
kappa.accept.rate = sum(kappa.count[(burn+1):nsim][seq(5,nsim-burn,by=5)])/length(kappa.count[(burn+1):nsim][seq(5,nsim-burn,by=5)]) #M-H acceptance rate
beta.accept.rate = sum(beta.count[(burn+1):nsim][seq(5,nsim-burn,by=5)])/length(beta.count[(burn+1):nsim][seq(5,nsim-burn,by=5)]) #M-H acceptance rate

beta.posterior.mean <- mean(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
beta.posterior.var <- var(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
beta.posterior.sd <- sqrt(beta.posterior.var)

kappa.posterior.mean <- mean(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
kappa.posterior.var <- var(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
kappa.posterior.sd <- sqrt(kappa.posterior.var)

alpha.posterior.mean <- mean(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
alpha.posterior.var <- var(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)])
alpha.posterior.sd <- sqrt(alpha.posterior.var)

beta.95th.cred.int <- quantile(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],probs=c(.025,.975))
alpha.95th.cred.int <- quantile(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],probs=c(.025,.975))
kappa.95th.cred.int <- quantile(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],probs=c(.025,.975))

latent1 <- mean(latent[(burn+1):nsim,1][seq(5,nsim-burn,by=5)]) 
latent2 <- mean(latent[(burn+1):nsim,2][seq(5,nsim-burn,by=5)])
latent3 <- mean(latent[(burn+1):nsim,3][seq(5,nsim-burn,by=5)])

result <- cbind(c(" ","posterior mean","posterior variance","posterior standard deviation","95th percentile lower credible bound","95th percentile upper credible bound"),
                c("alpha", alpha.posterior.mean, alpha.posterior.var, alpha.posterior.sd, alpha.95th.cred.int[[1]], alpha.95th.cred.int[[2]]),
                c("beta", beta.posterior.mean, beta.posterior.var, beta.posterior.sd, beta.95th.cred.int[[1]], beta.95th.cred.int[[2]]),
                c("kappa", kappa.posterior.mean, kappa.posterior.var, kappa.posterior.sd, kappa.95th.cred.int[[1]], kappa.95th.cred.int[[2]]),
                c("beta.accept.rate",beta.accept.rate,"kappa.accept.rate",kappa.accept.rate," "," "),
                c(latent1,latent2,latent3,"","",""))

postscript(file="data1.3parameterv7_MH_RW.eps",width = 5, height = 7, horizontal=F)
par(mfrow=c(2,3))
plot(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Trace Plots for ",alpha)),ylab=expression(alpha),type="l")
plot(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Trace Plots for ",beta)),ylab=expression(beta),type="l")
plot(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Trace Plots for ",kappa)),ylab=expression(kappa),type="l")

hist(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",alpha)),xlab=expression(alpha),prob=T)
lines(density(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))

hist(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",beta)),xlab=expression(beta),prob=T)
lines(density(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))

hist(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",kappa)),xlab=expression(kappa),prob=T)
lines(density(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))
dev.off()

write.table(cbind(beta,alpha2,alpha,kappa,rbind(rep("NA",m),latent)),file="Sampling_Scheme_1_v7_beta.RW_kappa.MH_Simulated_Alpha_Beta_kappa.csv",sep=",",
            col.names=c("beta","alpha^2","alpha","kappa","latent1","latent2","latent3"),row.names=F)
write.table(result,file="Sampling_Scheme_1_v7_beta.RW_kappa.MH_Simulated_Result.csv",sep=",",col.names=F,row.names=F)

# create data frame for Surv
Surv.dat <- data.frame( rbind( cbind(data, rep(2,length(data))),
                               cbind(cens.data, rep(1,length(cens.data))) ))
names(Surv.dat) <- c("Time","Status")
Surv.dat$SurvObj <- with(Surv.dat,Surv(Time, Status == 2))

## Kaplan-Meier estimator. The "log-log" confidence interval is preferred.
km <- survfit(SurvObj ~ 1, data = Surv.dat, conf.type = "log-log",
              type = "kaplan-meier")

rel.func <- function(a,b,k,t){
  1 - pnorm( 1/a * ( t^(1-k)/sqrt(b) - sqrt(b)/t^k  ) )
}

x = seq(from = min(c(data,cens.data)), to = 45, by = 0.01)
mle.Ft = rel.func(alpha.hat, beta.hat, kappa.hat, x)
gibbs.Ft = rel.func(alpha.posterior.mean, beta.posterior.mean, 
                    kappa.posterior.mean, x)

ps.filename = paste("data1.reliability.eps")
ps.plot.title = paste("Reliability Plot")
postscript(file=ps.filename,width = 3, height = 4, horizontal=F)
plot(km, conf.int=F, mark.time=F, ylab="R(t)", xlab="t", lty=1, col="black",
     main=ps.plot.title, ylim=c(0,1), xlim=c(min(c(data,cens.data)),45))
lines(x, mle.Ft, lty=4, lwd=2, col = "red")
lines(x, gibbs.Ft, lty=5, lwd=2, col = "blue")
legend("topright",inset=0.02, 
       c("Kaplan-Meier","MLE", "Gibbs"), 
       lty = c(1,4,5), 
       lwd = c(1,2,2),
       col = c("black","red", "blue"),
       cex=0.6)
dev.off()

#Histogram and Reliability Plot
hist.rel.plot.infile = paste("data1_hist_rel_plot.eps")
postscript(file=hist.rel.plot.infile ,width = 5, height = 5, horizontal=F)
par(mfrow=c(2,2))
hist(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",alpha)),xlab=expression(alpha),prob=T)
lines(density(alpha[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))

hist(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",beta)),xlab=expression(beta),prob=T)
lines(density(beta[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))

hist(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],main=expression(paste("Histogram for ",kappa)),xlab=expression(kappa),prob=T)
lines(density(kappa[(burn+2):(nsim+1)][seq(5,nsim-burn,by=5)],adjust=2))

plot(km, conf.int=F, mark.time=F, ylab="R(t)", xlab="t", lty=1, col="black",
     main="Reliability Plot", ylim=c(0,1), xlim=c(min(c(data,cens.data)),45))
lines(x, mle.Ft, lty=4, lwd=2, col = "red")
lines(x, gibbs.Ft, lty=5, lwd=2, col = "blue")

legend("topright", bty="n", cex=0.6,
       c("Kaplan-Meier","MLE", "Gibbs"), 
       lty = c(1,4,5), 
       lwd = c(1,2,2),
       col = c("black","red", "blue"))
par(mfrow=c(1,1))
dev.off()

