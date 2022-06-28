##---------------------------------------- ##
## Test 3PL model 
## Bambirra and Goncalves - vanilla (no mixture)
##----------------------------------------##
## lambda[i] - discrimination parameter
## beta[i]   - difficulty paramter
## delta[i]  - guessing parameter
## theta[j]  - individial ability

code3PL <- nimbleCode({

  for(j in 1:N) {
    for(i in 1:I) {
      y[j, i] ~ dbern(pi[j, i])
      pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[i, j]
      probit(linearReg[i, j]) <-  lambda[i]*eta[j] - gamma[i]
    }
    zi[j] ~ dcat(prob[1:2])
  }  
  
  for(i in 1:I) {
    lambda[i] ~ T(dnorm(1, var = 3), 0,)   
    gamma[i] ~ dnorm(0,  var = 3)
    delta[i] ~ dbeta(4, 12)
  } 

  
  prob[1:2] ~ ddirch(alpha[1:2])
  ## Mixture component parameter drawn from the base measure
  for(j in 1:N) {
    eta[j] ~ dnorm(mu[j], var = s2[j])  
    mu[j] <- muTilde[zi[j]]                 
    s2[j] <- s2Tilde[zi[j]]   
  }

  for(m in 1:2) {
    muTilde[m] ~ dnorm(0, var = s2_mu)
    s2Tilde[m] ~ dinvgamma(nu1, nu2)
  }


})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(gamma   = rnorm(constants$I, 0, 1),
              lambda = runif(constants$I, 0.5, 3), 
              delta  = rbeta(constants$I, 4, 12),
              nu1 = 2.01, nu2 = 1.01, s2_mu = 2,
              alpha = rep(1, 2))

monitors <- c("gamma", "lambda", "delta",
              "zi", "muTilde", "s2Tilde")



## Test 3PL model - Bambirra and Goncalves