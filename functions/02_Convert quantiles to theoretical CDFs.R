
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define functions for converting quantiles to the CDF 
###'       (for the theoretical quantiles)
###'
###' Date: 
###' - 2020-03-27: created
###' - 2022-02-14: updated -> rescaling with sigma_tau
###' 
###' Author: JoonHo Lee (`jlee296@ua.edu`)
###' 
###'

###'######################################################################
###' 
###'  Set quantiles (tail probabilities)
###' 
###'

vec_qt <- c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95)



###'######################################################################
###' 
###' Q2CDF_Gaussian()
###' 
###' - Convert quantiles to the "Gaussian" CDF
###' 
###' 

Q2CDF_Gaussian <- function(vec_qt, sigma_tau = 0.15){
  
  # Get quantile and CDF columnns
  Quantile <- vec_qt
  CDF <- qnorm(vec_qt)*sigma_tau # rescaling with sigma_tau
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_T()
###' 
###' - Convert quantiles to the "Student T" CDF
###' 
###' 

Q2CDF_T <- function(vec_qt, nu = 5,  sigma_tau = 0.15){
  
  Quantile <- vec_qt
  sd <- sqrt(nu/(nu - 2))  # scale with SD to have unit variance
  CDF <- qt(Quantile, nu)*(1/sd)*sigma_tau # rescale unit var with sigma_tau
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_Mixure()
###' 
###' - Convert quantiles to the CDF of "a mixture of two Gaussian distributions"
###' - Why theoretical? Because we simulate the large number (20,000)
###' 
###' 

Q2CDF_Mixture <- function(vec_qt, n_size = 20000, 
                          delta = 4, eps = 0.3, ups = 1, 
                          sigma_tau = 0.15){
  
  ### Get quantiles
  Quantile <- vec_qt 
  
  
  ### Calculate PDF (not CDF)
  ind <- runif(n_size) < (1 - eps)
  
  a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
  
  PDF <- ind*rnorm(n_size, -eps*delta/a, sqrt(1/a^2)) + 
    (1 - ind)*rnorm(n_size, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  PDF_scaled <- PDF*sigma_tau
  
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF_scaled, vec_qt))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_ALD()
###' 
###' - Convert quantiles to the CDF of asymmetric Laplace distribution
###' 
###' 

Q2CDF_ALD <- function(vec_qt, n_size = 20000, 
                      mean = 0, var = 1, p = 0.1, 
                      sigma_tau = 0.15){
  
  ### Get quantiles
  Quantile <- vec_qt 
  
  
  ###' Generate location, scale, and skewness parameters of ALD distribution
  ###' to get E(theta) = 0 and Var(theta) = 1
  scale <- sqrt((2*p^2*var)/(1 + p^4))
  location <- mean - ((scale*(1/p - p))/sqrt(2))
  
  
  ### Calculate PDF (not CDF)
  PDF <- LaplacesDemon::ralaplace(n_size, location, scale, p)
  
  PDF_scaled <- PDF*sigma_tau
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF_scaled, vec_qt))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}


###'######################################################################
###' 
###' Q2CDF_SkewN()
###' 
###' - Convert quantiles to the CDF of "Skewed Normal" distribution
###' 
###' 

Q2CDF_SkewN <- function(vec_qt, n_size = 20000, 
                        mean = 0, var = 1, slant = 5, 
                        sigma_tau = 0.15){
  
  ### Get quantiles
  Quantile <- vec_qt 
  
  
  ###' Generate location and scale parameters of skewed normal
  ###' to get E(Y) = 0 and Var(Y) = 1
  delta <- slant/sqrt(1 + slant^2)
  scale <- sqrt(var/(1-(2*delta^2/pi)))
  location <- mean - scale*sqrt(2/pi)*delta
  
  
  ### Calculate PDF (not CDF)
  PDF <- sn::rsn(n = n_size, xi = location, omega = scale, 
                 alpha = slant)[1:n_size]
  
  PDF_scaled <- PDF*sigma_tau
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF_scaled, vec_qt))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}


###'#######################################################################
###'
###' Generate true CDF distributions from quantiles: A wrapper
###' 
###' `quantile_to_CDF()`
###' 
###' - true_G [chr]: c("Gaussian", "T", "Skew", "ALD", "Bimodal", "Mixed")
###' - sigma_tau [dbl]: c(0.05, 0.10, 0.15, 0.20, 0.25)
###' - vec_qt [vec]: c(.05, .10, .25, .50, .75, .90, .95)
###' 
###' Because now the standard deviations of true distribution is not `1`
###' but `sigma_tau = {0.05, 0.10, 0.15, 0.20, 0.25}` 
###' 
###'

quantile_to_CDF <- function(true_G = "Gaussian", 
                            sigma_tau = 0.15, 
                            vec_qt = c(.05, .10, .25, .50, .75, .90, .95)){
  
  
  if (true_G == "Gaussian"){ 
    
    Q2CDF_Gaussian(vec_qt, sigma_tau)
    
  } else if (true_G == "T"){   
    
    Q2CDF_T(vec_qt, nu = 5,  sigma_tau)
    
  } else if (true_G == "Skew"){
    
    Q2CDF_SkewN(vec_qt, n_size = 20000, 
                mean = 0, var = 1, slant = 5, 
                sigma_tau)
    
  } else if (true_G == "ALD"){
    
    Q2CDF_ALD(vec_qt, n_size = 20000, 
              mean = 0, var = 1, p = 0.1, 
              sigma_tau)
    
  } else if (true_G == "Bimodal"){
    
    Q2CDF_Mixture(vec_qt, n_size = 20000, 
                  delta = 4, eps = 0.3, ups = 1, 
                  sigma_tau)
    
  } else if (true_G == "Mixed"){
    
    Q2CDF_Mixture(vec_qt, n_size = 20000, 
                  delta = 5, eps = 0.3, ups = 2, 
                  sigma_tau)
  }
}

# ### Test the function
# vec_qt <- c(.05, .10, .25, .50, .75, .90, .95)
# 
# quantile_to_CDF(true_G = "ALD", 
#                 sigma_tau = 0.15, 
#                 vec_qt = vec_qt)


