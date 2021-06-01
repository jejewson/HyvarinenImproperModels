
half_Gauaaian_nu_kappa <- function(kappa, m_0, s2_0){
  return(as.numeric(kappa >=0)*dnorm(1/kappa^2, m_0, sqrt(s2_0))/(1 - pnorm(0, m_0, sqrt(s2_0)))*2/kappa^3)
}

dkappa_pior_nu_IG <- function(kappa, a_0, b_0){
  return(2*exp(lgamma(a_0 + 0.5) - lgamma(a_0) - 0.5*log(b_0))*dgamma(kappa^2, shape = a_0 + 0.5, rate = b_0))
}

pkappa_pior_nu_IG <- function(kappa, a_0, b_0){
  return(integrate(function(k){dkappa_pior_nu_IG(k, a_0, b_0)}, lower = 0, upper = kappa)$value)
}

target_prior_density <- function(a_0, b_0, target){
  return(((pkappa_pior_nu_IG(target_upper, a_0, b_0) - pkappa_pior_nu_IG(target_lower, a_0, b_0)) - target)^2)
}


### Selecting a NLP for $\nu$

#We seek a minimally informative default non-local prior for 
#$\nu = \frac{1}{\kappa^2}$ in it's parametrisation of the Tukey's loss.

#$\nu = 0$ recovers the Gaussian model and therefore an Inverse-Gamma is 
#a standard NLP in this scenario. $\kappa$ was a nice parameterisation as 
#$\kappa$ was the number of standard deviations we consider cutting off. 
#I guess it would be reasonable to say that we want $\kappa\in\{3, 5\}$ with 
#probability 0.95. The CDF of the inverse-gamma distributions is 
#$F(x) = \frac{\Gamma(\alpha, \beta/x)}{\Gamma(x)}$ where the numerator 
#is the incomplete gamma function. So I am interested in 
#$\{\alpha, \beta: F(3; \alpha, \beta) - F(1; \alpha, \beta) = 0.95\}$.

target_upper <- 3
target_lower <- 1

target_prob <-  0.95

IG_prior_params <- optim(par = c(log(1), log(1)), fn = function(theta){target_prior_density(a_0 = exp(theta[1]), b_0 = exp(theta[2]), target = target_prob)})

a_0_nu_NLP_select <- exp(IG_prior_params$par[1])
b_0_nu_NLP_select <- exp(IG_prior_params$par[2])

