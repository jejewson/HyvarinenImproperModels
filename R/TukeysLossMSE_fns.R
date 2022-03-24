
grad_mu_log_score_tukey_varThresh <- function(x, mu, sigma2, c){
  ind <- as.numeric(abs(x - mu) <= (c*sqrt(sigma2)))
  return(ind*((x - mu)/sigma2 - 2*(x - mu)^3/(sigma2^2*c^2) + (x - mu)^5/(sigma2^3*c^4)))
}

grad_sigma2_log_score_tukey_varThresh <- function(x, mu, sigma2, c){
  ind <- as.numeric(abs(x - mu) <= (c*sqrt(sigma2)))
  return(-0.5*1/sigma2 + ind*((x - mu)^2/(2*sigma2^2) - (x - mu)^4/(sigma2^3*c^2) + (x - mu)^6/(2*sigma2^4*c^4)))
}

grad_mu_grad_mu_log_score_tukey_varThresh <- function(x, mu, sigma2, c){
  ind <- as.numeric(abs(x - mu) <= (c*sqrt(sigma2)))
  return(ind*(-1/sigma2 + 6*(x - mu)^2/(sigma2^2*c^2) - 5*(x - mu)^4/(sigma2^3*c^4)))
}

grad_mu_grad_sigma2_log_score_tukey_varThresh <- function(x, mu, sigma2, c){
  ind <- as.numeric(abs(x - mu) <= (c*sqrt(sigma2)))
  return(ind*(-(x - mu)/(sigma2^2) + 4*(x - mu)^3/(sigma2^3*c^2) - 3*(x - mu)^5/(sigma2^4*c^4)))
}

grad_sigma2_grad_sigma2_log_score_tukey_varThresh <- function(x, mu, sigma2, c){
  ind <- as.numeric(abs(x - mu) <= (c*sqrt(sigma2)))
  return(0.5*1/(sigma2^2) + ind*(- (x - mu)^2/(sigma2^3) + 3*(x - mu)^4/(sigma2^4*c^2) - 2*(x - mu)^6/(sigma2^5*c^4)))
}

J_log_score_tukey_varThresh <- function(g_x, mu, sigma2, kappa){
  Exp_grad_mu_grad_mu_log_score_tukey_varThresh <- integrate(function(x){grad_mu_grad_mu_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)*g_x(x)}, lower = -Inf, upper = Inf)
  Exp_grad_mu_grad_sigma2_log_score_tukey_varThresh <- integrate(function(x){grad_mu_grad_sigma2_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)*g_x(x)}, lower = -Inf, upper = Inf)
  Exp_grad_sigma2_grad_sigma2_log_score_tukey_varThresh <- integrate(function(x){grad_sigma2_grad_sigma2_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)*g_x(x)}, lower = -Inf, upper = Inf)
  return(matrix(c(Exp_grad_mu_grad_mu_log_score_tukey_varThresh$value, Exp_grad_mu_grad_sigma2_log_score_tukey_varThresh$value, Exp_grad_mu_grad_sigma2_log_score_tukey_varThresh$value, Exp_grad_sigma2_grad_sigma2_log_score_tukey_varThresh$value), nrow = 2, ncol = 2))
  
}


K_log_score_tukey_varThresh <- function(g_x, mu, sigma2, kappa){
  Exp_grad_mu_log_score_tukey_varThresh2 <- integrate(function(x){grad_mu_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)^2*g_x(x)}, lower = -Inf, upper = Inf)
  Exp_grad_mu_log_score_tukey_varThresh_grad_sigma2_log_score_tukey_varThresh <- integrate(function(x){grad_mu_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)*grad_sigma2_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)*g_x(x)}, lower = -Inf, upper = Inf)
  Exp_grad_sigma2_log_score_tukey_varThresh2 <- integrate(function(x){grad_sigma2_log_score_tukey_varThresh(x, mu, sigma2, c = kappa)^2*g_x(x)}, lower = -Inf, upper = Inf)
  return(matrix(c(Exp_grad_mu_log_score_tukey_varThresh2$value, Exp_grad_mu_log_score_tukey_varThresh_grad_sigma2_log_score_tukey_varThresh$value, Exp_grad_mu_log_score_tukey_varThresh_grad_sigma2_log_score_tukey_varThresh$value, Exp_grad_sigma2_log_score_tukey_varThresh2$value), nrow = 2, ncol = 2))
  
}


MSE_mu_log_score_tukey_varThresh <- function(n, g_x, kappa, mu_star = 0, theta_kappa_par){

  ## Sandwich matricies
  Jm1_eval <- solve(J_log_score_tukey_varThresh(g_x, mu = theta_kappa_par[1], sigma2 = theta_kappa_par[2], kappa))
  K_eval <- K_log_score_tukey_varThresh(g_x, mu = theta_kappa_par[1], sigma2 = theta_kappa_par[2], kappa)
  sandwich_cov <- Jm1_eval%*%K_eval%*%Jm1_eval
  
  return((theta_kappa_par[1] - mu_star)^2 + 1/n*sandwich_cov[1,1])
}


