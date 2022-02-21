
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define functions for performance evaluator (loss functions)
###'
###' Date: 2021-09-12
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###'

###'#######################################################################'
###'
###' `MSEL()`: the mean squared error loss 
###' 
###'  for the individual site-specific effect parameters tau_j's
###'
###'

MSEL <- function(RE.ests, RE.truth) {
  
  mean((RE.ests - RE.truth)^2)

}



###'#######################################################################'
###'
###' `MSELR()`: the mean squared error loss of the ranks 
###' 
###'  based on tau_j's
###'
###'

# NOTE: typo in paper after equation 6; equation 12 clarifies
#  - "rank" is from smallest to largest: smallest is rank 1, largest is rank J

MSELR <- function(RE.ests, RE.truth) {
  
  mean((rank(RE.ests) - rank(RE.truth))^2)

}



###'#######################################################################'
###'
###' `MSELP()`: the mean squared error loss of the percentiles 
###' 
###'  based on tau_j's
###'
###'

# NOTE: We have to rescale the MSELR estimates by 1/N^2 to get MSELP

MSELP <- function(RE.ests, RE.truth) {
  
  MSELR <- mean((rank(RE.ests) - rank(RE.truth))^2)
  J <- length(RE.ests)
  MSELR/(J^2)
  
}


###'#######################################################################'
###'
###' `ISEL()`: the integrated squared error loss
###' `IAEL()`: the integrated absolute error loss 
###'
###'  => the difference in the ECDFs across the support
###'     ECDF - Empirical Cumulative Distribution Function
###'
###'


ISEL <-  function( RE.ests, RE.truth ) {
  
  # Grid support of (RE.ests, RE.truth) into 400 pieces
  xs <- seq( min( RE.truth, RE.ests ), max( RE.truth, RE.ests ), length.out = 400 )
  delta <- xs[[2]] - xs[[1]]
  
  s1 <- ecdf( RE.ests )
  s2 <- ecdf( RE.truth )
  
  ISEL <- sum( (s1(xs) - s2(xs))^2 * delta )
  ISEL
}


ISEL_hhsim <- function( RE.ests, RE.truth ) {
  
  J <- length(RE.ests)
  
  ISEL_hhsim <- ensrisk(J, RE.truth, RE.ests, 0)$isel 
  ISEL_hhsim
}


IAEL = function( RE.ests, RE.truth ) {
  
  # Grid support of (RE.ests, RE.truth) into 400 pieces
  xs <- seq( min( RE.truth, RE.ests ), max( RE.truth, RE.ests ), length.out = 400 )
  delta <- xs[[2]] - xs[[1]]
  
  s1 <- ecdf( RE.ests )
  s2 <- ecdf( RE.truth )
  
  IAEL <- sum( abs(s1(xs) - s2(xs)) * delta )
  IAEL
}



###'#######################################################################'
###'
###' `KS_dist()`: Kolmogorov–Smirnov distance measures 
###' 
###'  for different forms of our density estimates
###'
###'

###' Looking at the KS distance of the empirical distribution vs. 
###' the true "empirical" distribution of the sample of sites we have.

KS_dist <- function( RE.ests, RE.truth ) {
  
  x <- sort( RE.ests )
  y <- sort( RE.truth )
  
  n.x <- length(x)
  n.y <- length(y)
  
  w <- c(x, y)
  
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  
  if (length(unique(w)) < (n.x + n.y)) {
    
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
    
  }
  
  KS_dist <- max( abs(z) )
  
  KS_dist
}



###'#######################################################################'
###'
###' Bias in skewness: 
###' 
###' skew(tau_j estimates) - skew(tau_j_true)
###' 
###' - `skew_true`
###' - `skew_est`
###' - `skewish_true`
###' - `skewish_est`
###' 
###'

# (1) The standard definition of skew:
skew <- function( X ) {
  
  mean( (X-mean(X))^3 ) / (sd(X)^3)

}

#' (2) a slightly more robust version:
#'     Roughly, skewish:skew::median:mean.
#'     NOTE: skewish = (n-1)/n if sign(X) is always +1
skewish <- function( X ) {
  
  mean( (X-mean(X))^2*sign(X) ) / (sd(X)^2)
  
}


