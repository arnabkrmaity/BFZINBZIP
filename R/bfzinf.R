#' This function computes log(Bayes factor) under Jefferey's prior
#' 
#' 
#' 
#' @param y a data vector
#'
#' @details NB: negative Binomial, ZINB: zero inflated negative Binomial,
#'          ZIP: zero inflated Poisson
#'          
#' @references Pramanik, P., and Maity, A. K. (2023) Bayes Factor of Zero Inflated Models under Jefferey's Prior.
#' 
#' @return 
#' \item{lbfNBP}{log(Bayes factor) of NB against Poisson}
#' \item{lbfZINBNB}{log(Bayes factor) of ZINB against NB} 
#' \item{lbfZIPP}{log(Bayes factor) of ZIP against Poisson}
#' \item{lbfZINBZIP}{log(Bayes factor) of ZINB against ZIP}
#' \item{propzero}{Porportion of Zeros present in the data}
#' 
#' 
#' @examples 
#' y <- rpois(n = 100, lambda = 5)
#' bfzinf(y)
#' 
#' @export
#' 


bfzinf <- function(y)
  # Returns the followings:
  # log(Bayes factor) of NB vs Poisson
  # log(Bayes factor) of ZINB vs NB
  # log(Bayes factor) of ZIP vs Poisson
  # log(Bayes factor) of ZINB vs ZIP
  # proportion of Zeros present in the data
  
{  
  n <- n_val <- length(y)
  
  ########################################################
  ######### Bayes factor of NB vs Poisson ################
  ########################################################
  
  gamma_val <- 1.001  # Empirically this works the best
  phi   <- sum(y)
  
  first  <- (phi + 0.5) * log(n) + lgamma(n * gamma_val) - 0.5 * log(gamma_val) - lgamma(phi + n * gamma_val + 0.5)
  second <- (n * gamma_val + 0.5) * sum(lgamma(y + gamma_val) - lgamma(gamma_val)) 
  third  <- (n * gamma_val - 0.5) * sum(- log(gamma(y + 1)))
  
  
  lbf1 <- first + second + third
  
  
  
  ########################################################
  ######### Bayes factor of ZINB vs NB ###################
  ########################################################
  
  
  tau_val <- 1.25  # Empirically this works the best
  # In Param's paper \tau is \gamma
  
  k_val = sum(y == 0)
  s_val = sum(y)
  m1_val = s_val - 0.5
  m2_val = s_val - 0.5
  n2_val = n_val*tau_val - 1
  
  nblr=array(0,c(n_val,1))
  zinblr=array(0,c(n_val,1))
  x2 = array(0, c(k_val+1, 1)) 
  
  lbeta2_val <- lbeta(m2_val + 1, n2_val + 1)
  lconstant_val <- lgamma(k_val + 1) - lgamma(n_val + 2)
  
  
  
  for (j_val in 1:(k_val+1)) 
  {
    n1_val = (n_val-(j_val-1))*tau_val - 1 
    lbeta1_val = lbeta(m1_val + 1, n1_val + 1)
    
    x2[j_val, ] <- (lgamma(n_val - (j_val - 1) + 1) - lgamma(k_val - (j_val - 1) + 1)) -
      (- j_val + 1) * tau_val * log(tau_val) + lbeta1_val
  }
  
  lbf2 <- matrixStats::logSumExp(x2) + lconstant_val - lbeta2_val
  
  
  
  
  ########################################################
  ######### Bayes factor of ZIP vs Poisson ###############
  ########################################################
  
  
  phi   <- sum(y)
  omega <- sum(y == 0)
  x <- array(0, c(omega + 1, 1)) 
  
  lconst <- lgamma(omega + 1) - lgamma(n + 2)
  
  for (j in 1:(omega + 1)) 
  {
    x[j, ] <- (lgamma(n - (j - 1) + 1)) - lgamma(omega - (j - 1) + 1) - 
      (phi + 0.5) * log(1 - (j - 1)/n)
  }
  
  lbf3 <- lconst + matrixStats::logSumExp(x)
  
  
  
  ########################################################
  ######### Bayes factor of ZINB vs ZIP ##################
  ########################################################
  
  
  zero.id <- which(y == 0)
  k_val = sum(y == 0)  # this is \omega in Param's paper
  zero.sample.id <- sample(zero.id, size = ceiling(k_val * 0.85), replace = FALSE, prob = NULL)
  y2 <- y[-zero.sample.id]
  y1 <- c(y2, y2)
  n1 <- length(y1)
  
  gamma_val <- 1.001  # Empirically this works the best
  phi   <- sum(y1)
  
  first  <- (phi + 0.5) * log(n1) + lgamma(n1 * gamma_val) - 0.5 * log(gamma_val) - lgamma(phi + n1 * gamma_val + 0.5)
  second <- (n1 * gamma_val + 0.5) * sum(lgamma(y1 + gamma_val) - lgamma(gamma_val)) 
  third  <- (n1 * gamma_val - 0.5) * sum(- log(gamma(y1 + 1)))
  
  lbf.test <- first + second + third
  
  
  tau_val <- ifelse(lbf.test >= log(3.2), yes = 1.25, no = 0.9)
  
  
  s_val = sum(y)       # this is \phi in Param's paper
  
  lfirst.const  <- -0.5 * log(tau_val)
  lsecond.const <- (n * tau_val + 0.5) * sum(lgamma(y + tau_val) - lgamma(tau_val))
  lthird.const  <- (n * tau_val - 0.5) * sum(-lgamma(y + 1))
  
  x <- array(0, c(k_val + 1, 1)) 
  
  for (j_val in 1:(k_val+1)) 
  {
    
    x[j_val, ] <- (lgamma((n - (j_val - 1)) * tau_val) - (s_val + 0.5) * log(n - (j_val - 1)) -
                     lgamma(s_val + (n - (j_val -1)) * tau_val + 0.5))
  }
  
  lbf4 <- lfirst.const + lsecond.const + lthird.const + matrixStats::logSumExp(x)
  
  
  result = list ("lbfNBP" = lbf1, "lbfZINBNB" = lbf2, 
                 "lbfZIPP" = lbf3, "lbfZINBZIP" = lbf4,
                 "propzero" = k_val/n)
  
  return(result)
}