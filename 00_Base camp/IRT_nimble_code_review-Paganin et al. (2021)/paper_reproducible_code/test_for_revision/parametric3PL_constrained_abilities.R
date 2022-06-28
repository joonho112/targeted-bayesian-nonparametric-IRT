##---------------------------------------- ##
## 3PL model - parametric
##----------------------------------------##
## lambda[i] - discrimination parameter
## beta[i]   - difficulty paramter
## delta[i]  - guessing parameter
## theta[j]  - individial ability

code3PL <- nimbleCode({
  for(i in 1:I) {
    for(j in 1:N) {
      y[j, i] ~ dbern(pi[j, i])
      pi[j, i] <- delta[i] + (1 - delta[i]) * linearReg[i, j]
      probit(linearReg[i, j]) <-  lambda[i]*eta[j] - gamma[i]
    }
  }  
  
  for(i in 1:I) {
    lambda[i] ~ T(dnorm(1, var = 3), 0,)   
    gamma[i] ~ dnorm(0,  var = 3)
    delta[i] ~ dbeta(4, 12)
  } 
  
  for(j in 1:N) {
    eta[j] ~ dnorm(0, 1)
  }

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(gamma   = rnorm(constants$I, 0, 1),
              lambda = runif(constants$I, 0.5, 3), 
              delta  = rbeta(constants$I, 4, 12)
              )

monitors <- c("gamma", "lambda", "delta")


## Test 3PL model - Babirra and Goncalves