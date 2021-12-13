library("DPpackage")

################################################
# simulated data
################################################
set.seed(0)

dtrue <- function(grid,x)
{
    	exp(-2*x)*dnorm(grid,mean=x,sd=sqrt(0.01))+
  (1-exp(-2*x))*dnorm(grid,mean=x^4,sd=sqrt(0.04))
} 

mtrue <- function(x)
{
       exp(-2*x)*x+(1-exp(-2*x))*x^4
} 

    nrec <- 500
    x <- runif(nrec)
    y1 <- x + rnorm(nrec, 0, sqrt(0.01))
    y2 <- x^4 + rnorm(nrec, 0, sqrt(0.04))
    u <- runif(nrec)
    prob <- exp(-2*x)
    y <- ifelse(u<prob,y1,y2)

################################################
# prior information - WDDP
################################################
w <- cbind(y,x)  
wbar <- apply(w,2,mean)
wcov <- var(w)

prior<-list(a0=10,
			    b0=1,
			    nu1=4,
			    nu2=4,
			s2=0.5*wcov,
			    m2=wbar,
			    psiinv2=2*solve(wcov),
			    tau1=6.01,
			    tau2=2.01)

################################################
# mcmc specification
################################################
mcmc <- list(nburn=5000,
			 nsave=5000,
			 nskip=3,
			 ndisplay=1000)

################################################
# covariate values where the density
# and mean function is evaluated
################################################
xpred <- seq(0,1,0.02)	

################################################
# fiiting the model - WDDP
################################################
fitWDDP <- DPcdensity(y=y,x=x,xpred=xpred,ngrid=100, 
  		  compute.band=TRUE,
  		  type.band="HPD",
  		  prior=prior, 
  		  mcmc=mcmc, 
  		  state=NULL, 
  		  status=TRUE)

################################################
# plots
################################################

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitWDDP$grid,fitWDDP$densp.h[6,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitWDDP$grid,fitWDDP$densp.l[6,],lwd=3,type="l",lty=2)
lines(fitWDDP$grid,fitWDDP$densp.m[6,],lwd=3,type="l",lty=1)
lines(fitWDDP$grid,dtrue(fitWDDP$grid,xpred[6]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitWDDP$grid,fitWDDP$densp.h[13,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitWDDP$grid,fitWDDP$densp.l[13,],lwd=3,type="l",lty=2)
lines(fitWDDP$grid,fitWDDP$densp.m[13,],lwd=3,type="l",lty=1)
lines(fitWDDP$grid,dtrue(fitWDDP$grid,xpred[13]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitWDDP$grid,fitWDDP$densp.h[25,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitWDDP$grid,fitWDDP$densp.l[25,],lwd=3,type="l",lty=2)
lines(fitWDDP$grid,fitWDDP$densp.m[25,],lwd=3,type="l",lty=1)
lines(fitWDDP$grid,dtrue(fitWDDP$grid,xpred[25]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitWDDP$grid,fitWDDP$densp.h[39,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitWDDP$grid,fitWDDP$densp.l[39,],lwd=3,type="l",lty=2)
lines(fitWDDP$grid,fitWDDP$densp.m[39,],lwd=3,type="l",lty=1)
lines(fitWDDP$grid,dtrue(fitWDDP$grid,xpred[39]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitWDDP$grid,fitWDDP$densp.h[45,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitWDDP$grid,fitWDDP$densp.l[45,],lwd=3,type="l",lty=2)
lines(fitWDDP$grid,fitWDDP$densp.m[45,],lwd=3,type="l",lty=1)
lines(fitWDDP$grid,dtrue(fitWDDP$grid,xpred[45]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(x,y,xlab="x",ylab="y",main="")
lines(xpred,fitWDDP$meanfp.m,type="l",lwd=3,lty=1)
lines(xpred,fitWDDP$meanfp.l,type="l",lwd=3,lty=2)
lines(xpred,fitWDDP$meanfp.h,type="l",lwd=3,lty=2)
    lines(xpred,mtrue(xpred),col="red",lwd=3)

################################################
# prior information - LDDP
################################################
library("splines")
W <- cbind(rep(1,nrec),bs(x,df=6,Boundary.knots=c(0,1)))
S0 <- 1000*solve(t(W)%*%W)
m0 <- solve(t(W)%*%W)%*%t(W)%*%y

prior<-list(a0=10,
			    b0=1,
			    m0=m0,
			    S0=S0,
			    tau1=6.01,
			    taus1=6.01,  
			    taus2=2.01,
			    nu=9,
			    psiinv=solve(S0))

################################################
# covariate values where the density
# and mean function is evaluated
################################################
xpred <- seq(0,1,0.02)	
Wpred <- cbind(rep(1,length(xpred)),predict(bs(x,df=6,Boundary.knots=c(0,1)),xpred))

################################################
# fiiting the model
################################################

fitLDDP <- LDDPdensity(formula=y~W-1,zpred=Wpred,
  		   ngrid=100, 
					   compute.band=TRUE,
  		   type.band="HPD",
  		   prior=prior, 
  		   mcmc=mcmc, 
  		   state=NULL, 
  		   status=TRUE)


################################################
# plots - LDDP
################################################

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitLDDP$grid,fitLDDP$densp.h[6,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitLDDP$grid,fitLDDP$densp.l[6,],lwd=3,type="l",lty=2)
lines(fitLDDP$grid,fitLDDP$densp.m[6,],lwd=3,type="l",lty=1)
lines(fitLDDP$grid,dtrue(fitLDDP$grid,xpred[6]),lwd=3,
      type="l",lty=1,col="red")


par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitLDDP$grid,fitLDDP$densp.h[13,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitLDDP$grid,fitLDDP$densp.l[13,],lwd=3,type="l",lty=2)
lines(fitLDDP$grid,fitLDDP$densp.m[13,],lwd=3,type="l",lty=1)
lines(fitLDDP$grid,dtrue(fitLDDP$grid,xpred[13]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitLDDP$grid,fitLDDP$densp.h[25,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitLDDP$grid,fitLDDP$densp.l[25,],lwd=3,type="l",lty=2)
lines(fitLDDP$grid,fitLDDP$densp.m[25,],lwd=3,type="l",lty=1)
lines(fitLDDP$grid,dtrue(fitLDDP$grid,xpred[25]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitLDDP$grid,fitLDDP$densp.h[39,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitLDDP$grid,fitLDDP$densp.l[39,],lwd=3,type="l",lty=2)
lines(fitLDDP$grid,fitLDDP$densp.m[39,],lwd=3,type="l",lty=1)
lines(fitLDDP$grid,dtrue(fitLDDP$grid,xpred[39]),lwd=3,
      type="l",lty=1,col="red")


par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(fitLDDP$grid,fitLDDP$densp.h[45,],lwd=3,type="l",lty=2,
     main="",xlab="y",ylab="f(y|x)",ylim=c(0,4))
lines(fitLDDP$grid,fitLDDP$densp.l[45,],lwd=3,type="l",lty=2)
lines(fitLDDP$grid,fitLDDP$densp.m[45,],lwd=3,type="l",lty=1)
lines(fitLDDP$grid,dtrue(fitLDDP$grid,xpred[45]),lwd=3,
      type="l",lty=1,col="red")

par(cex=1.5,mar=c(4.1, 4.1, 1, 1))
plot(x,y,xlab="x",ylab="y",main="")
lines(xpred,fitLDDP$meanfp.m,type="l",lwd=3,lty=1)
lines(xpred,fitLDDP$meanfp.l,type="l",lwd=3,lty=2)
lines(xpred,fitLDDP$meanfp.h,type="l",lwd=3,lty=2)
    lines(xpred,mtrue(xpred),col="red",lwd=3)


