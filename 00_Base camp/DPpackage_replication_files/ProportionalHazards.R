library("DPpackage")

################################################
# data
################################################
d <- read.table("kidney.txt") 

################################################
# function to make a row with '1' at ind
################################################
onv <- function(ind,len)
{
     onv <- rep(0,len)
     onv[ind] <- 1
     return(onv)
}

################################################
# Create data to fit Cox model using 
# Poisson likelihood for piecewise 
# exponential model.
################################################
ewdat <- matrix(1:(38*2*2),nrow=38*2,ncol=2)
t <- rep(0,38*2)
elta <- tt
or(i in 1:38)
 
   newdat[i*2-1,1] <- d[i,1]
   newdat[i*2-1,2] <- d[i,6]
   newdat[i*2  ,1] <- d[i,1]
   newdat[i*2  ,2] <- d[i,6]
   tt[i*2-1] <- d[i,2]
   delta[i*2-1] <- d[i,3]
   tt[i*2] <- d[i,4]
   delta[i*2] <- d[i,5]


 <- NULL
at <- NULL
ot <- 0
 <- ncol(newdat)
ff <- NULL
 <- length(tt)
ntervals <- 10
utpoint <- quantile(tt,(1:intervals)/intervals,names=FALSE)

or(i in 1:n)

   tot <- tot+1
   mat <- matrix(append(mat,c(newdat[i,1:p],onv(1,intervals))),
		       c(p+intervals,tot))
   off <- append(off,min(cutpoint[1],tt[i]))
   if(tt[i]<=cutpoint[1] && delta[i]==1)
   {
      y <- append(y,1)
   }
   else 
   {
      y <- append(y,0)
   } 
   for(j in 1:(intervals-1)) 
   {
       if(tt[i]>cutpoint[j]) 
       {
	  off <- append(off,min(cutpoint[j+1],
				tt[i])-cutpoint[j])
	  tot <- tot+1
	  mat <- matrix(append(mat,c(newdat[i,1:p],
			      onv(j+1,intervals))),
			      c(p+intervals,tot))
	  if(tt[i] <= cutpoint[j+1] && delta[i]==1)
	  {
	     y <- append(y,1)
	  } 
	  else
	  {
	     y <- append(y,0)
	  }
       }
   }
}
mat <- t(mat)
id <- mat[,1]
gender <- mat[,2]
loghazard <- mat[,3:12]

################################################
# PQL estimation
################################################
library("MASS")
fit0 <- glmmPQL(fixed=y~gender+loghazard-1+offset(log(off)),
        			random=~1|id,family=poisson(log))

################################################
# prior
################################################
beta0 <- fit0$coefficients$fixed
Sbeta0 <- vcov(fit0) 

prior <- list(M=5,
	      a0=1,
	      b0=1,
	      nu0=3,
	      tinv=diag(1,1),
	      mu=rep(0,1),
	      beta0=beta0,
	      Sbeta0=Sbeta0,
	      frstlprob=TRUE)

################################################
# starting values from PQL estimation
################################################
beta <- fit0$coefficients$fixed
b <- as.vector(fit0$coefficients$random$id)
mu <- rep(0,1)
sigma <- getVarCov(fit0)[1,1]

state <- list(alpha=1,
	      beta=beta,
	      b=b,
	      mu=mu,
	      sigma=sigma)

################################################
# mcmc
################################################
mcmc <- list(nburn=5000,
        		 nsave=5000,
	     nskip=19,
	     ndisplay=1000,
	     tune3=1.5)

################################################
# fitting the model
################################################
fitPT <- PTglmm(fixed=y~gender+loghazard,
		offset=log(off),
		random=~1|id,
		family=poisson(log),  
		prior=prior,
		mcmc=mcmc,
		state=state,
		status=FALSE)

################################################
# posterior inferences
################################################
summary(fitPT,hpd=FALSE)

################################################
# frailties density estimate
################################################
predPT <- PTrandom(fitPT,predictive=TRUE,gridl=c(-2.3,2.3))

#########################################################
# plots and output
#########################################################

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(predPT,ask=FALSE,lwd=3)

# survival curve function
# inputs: time, covariates, vector w/ (beta, log-hazards, predictive frailty), # intervals
Sc <- function(tt,x,beta.h.f,intervals)
{
	p <- length(x)
	j <- 1
	s <- min(cutpoint[1],tt)*exp(beta.h.f[p+1])
	while (tt>cutpoint[j] && j<intervals)     
	{
	   s <- s+(min(cutpoint[j+1],tt)-cutpoint[j])*exp(beta.h.f[p+j])
	   j=j+1
	}
	exp(-s*exp((x%*%beta.h.f[1:p])[1,1]+beta.h.f[p+intervals+1]))
}

# combine beta, log-hazard values on partition, and predictive (39th) frailty
all <- cbind(fitPT$save.state$thetasave[,2:(p+intervals)],fitPT$save.state$randsave[,39])

# men have covariate vector c(0)
time <- seq(0,300,5)
u <- time
m <- time
l <- time
dummy <- rep(0,5000)
for (i in 1:length(m)) 
{ 
	for(j in 1:5000) 
	{ 
	    dummy[j] <- Sc(time[i],c(0),all[j,],10) 
	}
    q <- quantile(dummy, c(0.025,0.5,0.975),names=FALSE)
	l[i] <- q[1]
	m[i] <- q[2]
	u[i] <- q[3]
}

plot(time,m,type="l",ylab="survival",ylim=c(0,1),lwd=3) 
lines(time,l,lty=2,lwd=3) 
lines(time,u,lty=2,lwd=3)

# women have covariate vector c(1)
time <- seq(0,300,5)
u <- time
m <- time
l <- time
dummy <- rep(0,5000)
for (i in 1:length(m)) 
{ 
	for(j in 1:5000) 
	{ 
  	    dummy[j] <- Sc(time[i],c(1),all[j,],10) 
	}
	q <- quantile(dummy, c(0.025,0.5,0.975),names=FALSE)
	l[i] <- q[1]
	m[i] <- q[2]
	u[i] <- q[3]
}

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(time,m,type="l",ylab="survival",ylim=c(0,1),lwd=3) 
lines(time,l,lty=2,lwd=3) 
lines(time,u,lty=2,lwd=3)


