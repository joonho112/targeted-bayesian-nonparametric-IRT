library("DPpackage")

################################################
# data
################################################

dat <- read.csv(file="base500.csv",header=TRUE) 
attach(dat)

y <- cbind(bin1,bin2,bin3,bin4,bin5,bin6,bin7,bin8,bin9,bin10,
           bin11,bin12,bin13,bin14,bin15,bin16,bin17,bin18,bin19,bin20,
           bin21,bin22,bin23,bin24,bin25,bin26,bin27,bin28,bin29,bin30,
           bin31,bin32,bin33,bin34,bin35,bin36,bin37,bin38,bin39,bin40,
           bin41,bin42,bin43,bin44,bin45)
  
types <- as.factor(DDCIA06)
gender <- GENERO

################################################
# prediction's design matrix 
# columns: - intercept
#	   - 3 dummy type of school
#	   - gender indicator 
################################################
zpred <- matrix(c(1,0,0,0,0,
        	  1,1,0,0,0,
        	  1,0,1,0,0,
        	  1,0,0,1,0,
        	  1,0,0,0,1,
        	  1,1,0,0,1,
        	  1,0,1,0,1,
        	  1,0,0,1,1),nrow=8,ncol=5,byrow=TRUE)

################################################
# prior information
################################################
prior <- list(alpha=1, 
             beta0=rep(0,44),
             Sbeta0=diag(1000,44),
             mu0=rep(0,5),
             S0=diag(100,5),
             tau1=6.01,
             taus1=6.01,
             taus2=2.01,
             nu=8,
             psiinv=diag(1,5))

################################################
# mcmc
################################################
mcmc <- list(nburn=5000,
			     nskip=3,
             ndisplay=1000,
             nsave=5000)

################################################
# fitting the model
################################################
fitLDDP <-  LDDPrasch(formula=y ~ types+gender,
        	      prior=prior,
        	      mcmc=mcmc,
        	      state=NULL,
        				  status=TRUE,
        	      zpred=zpred,
        	      grid=seq(-3,8,len=100),
        	      compute.band=TRUE)

###########################
# plots
###########################

# Boys - Public I
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[2,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[2,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[2,],lty=2,lwd=3,col=1)

# Girls - Public I
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[6,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[6,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[6,],lty=2,lwd=3,col=1)

# Boys - Public II
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[1,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[1,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[1,],lty=2,lwd=3,col=1)

# Girls - Public II
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[5,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[5,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[5,],lty=2,lwd=3,col=1)

# Boys - Private I
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[4,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[4,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[4,],lty=2,lwd=3,col=1)

# Girls - Private I
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[8,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[8,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[8,],lty=2,lwd=3,col=1)

# Boys - Private II
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[3,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[3,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[3,],lty=2,lwd=3,col=1)

# Girls - Private II
  par(cex=2,mar=c(4.1, 4.1, 1, 1))
  plot(fitLDDP$grid,fitLDDP$dens.m[7,],xlim=c(-2,7.3),ylim=c(0,0.6),type="l",lty=1,lwd=3,xlab=expression(theta),ylab="density",col=1)
  lines(fitLDDP$grid,fitLDDP$dens.u[7,],lty=2,lwd=3,col=1)
  lines(fitLDDP$grid,fitLDDP$dens.l[7,],lty=2,lwd=3,col=1)

