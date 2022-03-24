ISE_KDE_Gaussian <- function(z, h, bar_y = 0, S_y = 1){
  ISE <- integrate(f = function(x){(f_KDE_Gaussian_vect(x, data = z, sigma2 = h, w = 1) - dnorm(x, bar_y, sqrt(S_y)))^2}, lower = min(z), upper = max(z))
  return(ISE$value)
}

ISE_Gaussian_Gaussian <- function(bar_z, S_z, bar_y = 0, S_y = 1, min_z, max_z){
  ISE <- integrate(f = function(x){(dnorm(x, bar_z, sqrt(S_z)) - dnorm(x, bar_y, sqrt(S_y)))^2}, lower = min_z, upper = max_z)
  return(ISE$value)
}

Exp_ISE_Gaussian_Gaussian <- function(n, bar_y = 0, S_y = 1, N_MC){
  ISE_Gaussian_Gaussian_eval <- rep(NA, N_MC)
  z_sim <- array(NA, dim = c(n, N_MC))
  for(i in 1:N_MC){
    z_sim[, i] <- rnorm(n, bar_y, sqrt(S_y))
    ISE_Gaussian_Gaussian_eval[i] <- ISE_Gaussian_Gaussian(mean(z_sim[,i]), sqrt(var(z_sim[,i])), bar_y, S_y, min_z = -Inf, max_z = Inf)
  }
  return(ISE_Gaussian_Gaussian_eval)
}

Exp_ISE_KDE_Gaussian_h <- function(h, n, bar_y = 0, S_y = 1, N_MC){
  ISE_KDE_Gaussian_eval <- rep(NA, N_MC)
  z_sim <- array(NA, dim = c(n, N_MC))
  for(i in 1:N_MC){
    z_sim[, i] <- rnorm(n, bar_y, sqrt(S_y))
    ISE_KDE_Gaussian_eval[i] <- ISE_KDE_Gaussian(z_sim[,i], h, bar_y, S_y)
  }
  return(ISE_KDE_Gaussian_eval)
}

invgamma_select <- function(lower, upper, prob){
  optim_fn <- function(log_p){abs(prob - (pinvgamma(upper, shape = exp(log_p[1]), scale = exp(log_p[2])) - pinvgamma(lower, shape = exp(log_p[1]), scale = exp(log_p[2]))))}
  optim_par <- optim(par = c(log(1), log(1)), fn = optim_fn)
  return(exp(optim_par$par))
}

ISE_KDE_w_Gaussian <- function(z, h, w, bar_y = 0, S_y = 1){
  KDE_predictive_normaliser <- integrate(f = function(x){f_KDE_Gaussian_vect(x, data = z, sigma2 = h, w)}, lower = -Inf, upper = Inf)
  
  ISE <- integrate(f = function(x){(f_KDE_Gaussian_vect(x, data = z, sigma2 = h, w)/KDE_predictive_normaliser$value - dnorm(x, bar_y, sqrt(S_y)))^2}, lower = min(z), upper = max(z))
  return(ISE$value)
}

Exp_ISE_KDE_w_Gaussian_h <- function(h, w, n, bar_y = 0, S_y = 1, N_MC){
  ISE_KDE_Gaussian_eval <- rep(NA, N_MC)
  z_sim <- array(NA, dim = c(n, N_MC))
  for(i in 1:N_MC){
    z_sim[, i] <- rnorm(n, bar_y, sqrt(S_y))
    ISE_KDE_Gaussian_eval[i] <- ISE_KDE_w_Gaussian(z = z_sim[,i], h = h, w = w, bar_y, S_y)
  }
  return(ISE_KDE_Gaussian_eval)
}

Exp_ISE_KDE_Gaussian_hIG <- function(a_0, b_0, n, bar_y = 0, S_y = 1, N_MC){
  ISE_KDE_Gaussian_eval <- rep(NA, N_MC)
  h_sim <- rep(NA, N_MC)
  z_sim <- array(NA, dim = c(n, N_MC))
  for(i in 1:N_MC){
    z_sim[, i] <- rnorm(n, bar_y, sqrt(S_y))
    h_sim[i] <- rinvgamma(1, shape = a_0, scale = b_0)
    ISE_KDE_Gaussian_eval[i] <- ISE_KDE_Gaussian(z_sim[,i], h_sim[i], bar_y, S_y)
  }
  return(ISE_KDE_Gaussian_eval)
}

Exp_ISE_KDE_w_Gaussian_hIGExp <- function(a_0, b_0, lambda, n, bar_y = 0, S_y = 1, N_MC){
  ISE_KDE_Gaussian_eval <- rep(NA, N_MC)
  h_sim <- rep(NA, N_MC)
  w_sim <- rep(NA, N_MC)
  z_sim <- array(NA, dim = c(n, N_MC))
  for(i in 1:N_MC){
    z_sim[, i] <- rnorm(n, bar_y, sqrt(S_y))
    h_sim[i] <- rinvgamma(1, shape = a_0, scale = b_0)
    w_sim[i] <- rexp(1, lambda) # mean = 1/lambda
    ISE_KDE_Gaussian_eval[i] <- ISE_KDE_w_Gaussian(z = z_sim[,i], h = h_sim[i], w = w_sim[i], bar_y, S_y)
  }
  return(ISE_KDE_Gaussian_eval)
}

invgamma_exp_select <- function(n, N_MC){
  optim_fn <- function(log_p){mean(Exp_ISE_KDE_w_Gaussian_hIGExp(a_0 = 2, b_0 = exp(log_p[1]), lambda = exp(log_p[2]), n, bar_y = 0, S_y = 1, N_MC))}
  optim_par <- optim(par = c(log(invgamma_prior_set[2]), log(1)), fn = optim_fn)
  return(list("par" = exp(optim_par$par), "value" = optim_par$value))
}

