
### Hyvarinen score functions

#### log-score Gaussian (squared error-loss) ####

log_score_norm<- function(x, mu, sigma2,w){
  return(w*dnorm(x, mu, sqrt(sigma2),log=TRUE))
}

grad_log_score_norm <- function(x, mu, sigma2,w){
  return(-w*(x-mu)/sigma2)
}

Laplacian_log_score_norm <- function(x, mu, sigma2,w){
  return(-w*1/sigma2)
}

H_score_norm <- function(x, mu, sigma2, w){
  return(2*Laplacian_log_score_norm(x, mu, sigma2, w) + (grad_log_score_norm(x, mu, sigma2, w))^2)
}

#### Tukey's loss normal ####

tukey_loss <- function(x,c){
  ind <- as.numeric(abs(x) <= (c))
  return(ind*(x^2/2 - x^4/(2*c^2) + x^6/(6*c^4))+(1-ind)*(c^2/6))
}

log_score_tukey_varThresh <- function(x, mu, sigma2,c)
{
  return(-1/2*log(2*pi*sigma2)-(tukey_loss((x-mu)/sqrt(sigma2),c)))
}

### We need these functions without the indicators for the derivatives and Hessians given the abs approx

grad_log_score_tukey_varThresh_noind <- function(x, mu, sigma2,c)
{
  ind <- 1
  return(-(ind*((x-mu)/sigma2-2*(x-mu)^3/(c^2*sigma2^2)+(x-mu)^5/(c^4*sigma2^3))))
}

grad_log_score_tukey_varThresh <- function(x, mu, sigma2,c)
{
  ind <- as.numeric(abs(x-mu) <= (c*sqrt(sigma2)))
  return(ind * grad_log_score_tukey_varThresh_noind(x, mu, sigma2, c))
}

Laplacian_log_score_tukey_varThresh_noind <- function(x, mu, sigma2, c)
{
  ind <- 1
  return(-(ind*(1/sigma2-6*(x-mu)^2/(c^2*sigma2^2)+5*(x-mu)^4/(c^4*sigma2^3))))
}

Laplacian_log_score_tukey_varThresh <- function(x, mu, sigma2,c)
{
  ind <- as.numeric(abs(x-mu) <= (c*sqrt(sigma2)))
  return(ind * Laplacian_log_score_tukey_varThresh_noind(x, mu, sigma2, c))
  
}

H_score_tukey_varThresh_noind <- function(x, mu, sigma2, c)
{
  return(2*Laplacian_log_score_tukey_varThresh_noind(x, mu, sigma2,c) + (grad_log_score_tukey_varThresh_noind(x, mu, sigma2,c))^2)
}

H_score_tukey_varThresh <- function(x, mu, sigma2, c)
{
  return(2*Laplacian_log_score_tukey_varThresh(x, mu, sigma2,c) + (grad_log_score_tukey_varThresh(x, mu, sigma2,c))^2)
}

sigmoid <- function(x, k){
  #return(1 / (1 + exp( - k * x)))
  return(exp(- log1pexp(- k * x)))
}

abs_approx <- function(x, k){
  return(sqrt(x^2 + 1/k))
}

approx_indicator_abs_less3 <- function(x, c, k_abs, k_sigmoid){
  #return(sigmoid( - (sqrt((x)^2 + 1/(k_abs)) - c), k_sigmoid))
  return(sigmoid( - (abs_approx(x, k_abs) - c), k_sigmoid))
}


tukey_loss_absapprox <- function(x,c, k_abs, k_sigmoid){
  ind <- approx_indicator_abs_less3(x, c, k_abs, k_sigmoid);
  return(ind*(x^2/2 - x^4/(2*c^2) + x^6/(6*c^4))+(1-ind)*(c^2/6))
}

log_score_tukey_varThresh_absapprox <- function(x, mu, sigma2, c, k_abs, k_sigmoid)
{
  return(-1/2*log(2*pi*sigma2)-(tukey_loss_absapprox((x-mu)/sqrt(sigma2), c, k_abs, k_sigmoid)))
}

grad_log_score_tukey_varThresh_absapprox <- function(x, mu, sigma2, c, k_abs, k_sigmoid)
{
  ind <- approx_indicator_abs_less3((x - mu)/sqrt(sigma2), c, k_abs, k_sigmoid);
  return(-(ind*((x-mu)/sigma2-2*(x-mu)^3/(c^2*sigma2^2)+(x-mu)^5/(c^4*sigma2^3))))
  
}

Laplacian_log_score_tukey_varThresh_absapprox <- function(x, mu, sigma2, c, k_abs, k_sigmoid)
{
  ind <- approx_indicator_abs_less3((x - mu)/sqrt(sigma2), c, k_abs, k_sigmoid);
  return(-(ind*(1/sigma2-6*(x-mu)^2/(c^2*sigma2^2)+5*(x-mu)^4/(c^4*sigma2^3))))
  
}

H_score_tukey_varThresh_absapprox <- function(x, mu, sigma2, c, k_abs, k_sigmoid)
{
  return(2*Laplacian_log_score_tukey_varThresh_absapprox(x, mu, sigma2, c, k_abs, k_sigmoid) + (grad_log_score_tukey_varThresh_absapprox(x, mu, sigma2,c, k_abs, k_sigmoid))^2)
}

### Laplace Approximation Functions

log_laplace_approximation_marg_lik <- function(log_P_star, A, p){
  return(log_P_star + p/2*log(2*pi) - 1/2*log(det(A)))
}


### NIGG-prior Hessians

## Assumes independence # m is for minus
NIG_mlog_prior_regression <- function(beta, sigma2, a_0, b_0, beta_0, kappa_0){
  return(- dinvgamma(sigma2, shape = a_0, scale = b_0, log = TRUE) - sum(dnorm(beta, beta_0, sqrt(sigma2/kappa_0), log = TRUE)))
}

NIG_mlog_prior_Hessian_regression  <- function(beta, sigma2, a_0, b_0, beta_0, kappa_0){
  p <- length(beta)
  Hessian <- matrix(0, nrow = p + 1, ncol = p + 1)
  Hessian[1:p, 1:p] <- diag(kappa_0/sigma2, p)
  Hessian[p + 1, 1:p] <- - kappa_0*(beta - beta_0)/(sigma2^2)
  Hessian[1:p, p + 1] <- Hessian[p + 1, 1:p]
  Hessian[p + 1, p + 1] <- (2*b_0 + sum(kappa_0*(beta - beta_0)^2))/(sigma2^3) - (a_0 + (2 + p)/2)/(sigma2^2)
  
  return(Hessian)
}

inverse_gamma_mlog_prior_Hessian  <- function(kappa, a_0, b_0){
  return(- (a_0 + 1)/(kappa^2) + 2*b_0/(kappa^3))
}

NIGIG_kappa_mlog_prior_Hessian_regression <- function(beta, sigma2, kappa, a_0, b_0, beta_0, kappa_0, kappa_a_0, kappa_b_0){
  p <- length(beta)
  Hessian <- matrix(0, nrow = p + 2, ncol = p + 2)
  Hessian[1:(p+1), 1:(p+1)] <- NIG_mlog_prior_Hessian_regression(beta, sigma2, a_0, b_0, beta_0, kappa_0)
  Hessian[p + 2, p + 2] <- inverse_gamma_mlog_prior_Hessian(kappa, kappa_a_0, kappa_b_0)
  
  return(Hessian)
}

gamma_mlog_prior_Hessian  <- function(nu, a_0, b_0){
  return(-exp(-2*log(nu) + log(-(a_0 - 1)/(2) + b_0*sqrt(nu)/(4))))
}

NIGG_nu_mlog_prior_Hessian_regression <- function(beta, sigma2, nu, a_0, b_0, beta_0, kappa_0, nu_a_0, nu_b_0){
  p <- length(beta)
  Hessian <- matrix(0, nrow = p + 2, ncol = p + 2)
  Hessian[1:(p+1), 1:(p+1)] <- NIG_log_prior_Hessian_regression(beta, sigma2, a_0, b_0, beta_0, kappa_0)
  
  return(Hessian)
}

### Gaussian model Hessians

hess_H_score_norm_regression <- function(y, X, beta, sigma2){
  
  n <- nrow(X)
  p <- length(beta)
  
  Hess <- array(NA, dim = c(p + 1, p + 1, n))
  
  for(j in 1:p){
    for(k in 1:p){
      Hess[j, k, ] <- 2*X[,j]*X[,k]/(sigma2^2)
    }
    Hess[j, p + 1, ] <- 4*X[,j]*(y - X%*%beta)/(sigma2^3)
    Hess[p + 1, j, ] <- Hess[j, p + 1, ]
  }
  Hess[p + 1, p + 1,] <- -4/(sigma2^3) + 6*(y - X%*%beta)^2/(sigma2^4)
  
  return(Hess)
}

### Tukey's loss Hessians

# Reparametrising nu = 1/kappa^2 #

grad_H_score_Tukeys_regression_repar_noind <- function(y, X, beta, sigma2, nu){
  
  n <- nrow(X)
  p <- length(beta)
  
  ind <- 1
  
  kappa <- 1/sqrt(nu)
  
  grad_beta <- matrix(NA, nrow = p, ncol = n)
  for(j in 1:p){
    grad_beta[j,] <- ind*(-2*(12*X[,j]*(y - X%*%beta)/(kappa^2*sigma2^2) - 20*X[,j]*(y - X%*%beta)^3/(kappa^4*sigma2^3)) + 
                            (- 2*X[,j]*(y - X%*%beta)/(sigma2^2) + 16*X[,j]*(y - X%*%beta)^3/(kappa^2*sigma2^3) - 36*X[,j]*(y - X%*%beta)^5/(kappa^4*sigma2^4) + 32*X[,j]*(y - X%*%beta)^7/(kappa^6*sigma2^5) - 10*X[,j]*(y - X%*%beta)^9/(kappa^8*sigma2^6)))
  }
  
  grad_sigma2 <- ind*(-2*(- 1/sigma2^2 + 12*(y - X%*%beta)^2/(kappa^2*sigma2^3) - 15*(y - X%*%beta)^4/(kappa^4*sigma2^4)) + 
                        (- 2*(y - X%*%beta)^2/(sigma2^3) + 12*(y - X%*%beta)^4/(kappa^2*sigma2^4) - 24*(y - X%*%beta)^6/(kappa^4*sigma2^5) + 20*(y - X%*%beta)^8/(kappa^6*sigma2^6) - 6*(y - X%*%beta)^10/(kappa^8*sigma2^7)))
  
  grad_nu <- ind*(-2*(- 6*(y - X%*%beta)^2/(sigma2^2) + 10*nu*(y - X%*%beta)^4/(sigma2^3)) + 
                    (- 4*(y - X%*%beta)^4/(sigma2^3) + 12*nu*(y - X%*%beta)^6/(sigma2^4) - 12*nu^2*(y - X%*%beta)^8/(sigma2^5) + 4*nu^3*(y - X%*%beta)^10/(sigma2^6)))
  
  grad <- array(NA, dim = c(p + 2, n))
  
  grad[1:p, ] <- grad_beta
  grad[p + 1, ] <- grad_sigma2
  grad[p + 2, ] <- grad_nu
  
  return(grad)
}

grad_H_score_Tukeys_regression_repar <- function(y, X, beta, sigma2, nu){
  
  n <- length(y)
  p <- length(beta)
  
  ind <- array(NA, dim = c(p + 2, n))
  for(i in 1:(p + 2)){
    ind[i,] <- as.numeric((abs(y - X%*%beta) <= (sqrt(sigma2/nu))))
  }
  return(ind * grad_H_score_Tukeys_regression_repar_noind(y, X, beta, sigma2, nu))
}

hess_H_score_Tukeys_regression_repar_noind  <- function(y, X, beta, sigma2, nu){
  
  n <- nrow(X)
  p <- length(beta)
  
  ind <- 1
  
  kappa <- 1/sqrt(nu)
  
  Hess <- array(NA, dim = c(p + 2, p + 2, n))
  for(j in 1:p){
    for(k in 1:p){
      Hess[j,k,] <- ind*(-2*(-12*X[,j]*X[,k]/(kappa^2*sigma2^2) + 60*X[,j]*X[,k]*(y - X%*%beta)^2/(kappa^4*sigma2^3)) + 
                           (2*X[,j]*X[,k]/(sigma2^2) - 48*X[,j]*X[,k]*(y - X%*%beta)^2/(kappa^2*sigma2^3) + 180*X[,j]*X[,k]*(y - X%*%beta)^4/(kappa^4*sigma2^4) - 224*X[,j]*X[,k]*(y - X%*%beta)^6/(kappa^6*sigma2^5) + 90*X[,j]*X[,k]*(y - X%*%beta)^8/(kappa^8*sigma2^6)))
    }
    Hess[j, p + 1,] <- ind*(-2*(-24*X[,j]*(y - X%*%beta)/(kappa^2*sigma2^3) + 60*X[,j]*(y - X%*%beta)^3/(kappa^4*sigma2^4)) + 
                              (4*X[,j]*(y - X%*%beta)/(sigma2^3) - 48*X[,j]*(y - X%*%beta)^3/(kappa^2*sigma2^4) + 144*X[,j]*(y - X%*%beta)^5/(kappa^4*sigma2^5) - 160*X[,j]*(y - X%*%beta)^7/(kappa^6*sigma2^6) + 60*X[,j]*(y - X%*%beta)^9/(kappa^8*sigma2^7)))
    Hess[p + 1, j,] <- Hess[j, p + 1,]
    
    Hess[j, p + 2,] <- ind*(-2*(+12*X[,j]*(y - X%*%beta)/(sigma2^2) - 40*nu*X[,j]*(y - X%*%beta)^3/(sigma2^3)) + 
                              (16*X[,j]*(y - X%*%beta)^3/(sigma2^3) - 72*nu*X[,j]*(y - X%*%beta)^5/(sigma2^4) + 96*nu^2*X[,j]*(y - X%*%beta)^7/(sigma2^5) - 40*nu^3*X[,j]*(y - X%*%beta)^9/(sigma2^6)))
    Hess[p + 2, j,] <- Hess[j, p + 2,]
  }
  
  Hess[p + 1, p + 1,] <- ind*(-2*(2/sigma2^3 - 36*(y - X%*%beta)^2/(kappa^2*sigma2^4) + 60*(y - X%*%beta)^4/(kappa^4*sigma2^5)) + 
                                (6*(y - X%*%beta)^2/(sigma2^4) - 48*(y - X%*%beta)^4/(kappa^2*sigma2^5) + 120*(y - X%*%beta)^6/(kappa^4*sigma2^6) - 120*(y - X%*%beta)^8/(kappa^6*sigma2^7) + 42*(y - X%*%beta)^10/(kappa^8*sigma2^8)))
  
  Hess[p + 1, p + 2,] <- ind*(-2*(12*(y - X%*%beta)^2/(sigma2^3) - 30*nu*(y - X%*%beta)^4/(sigma2^4)) + 
                                (12*(y - X%*%beta)^4/(sigma2^4) - 48*nu*(y - X%*%beta)^6/(sigma2^5) + 60*nu^2*(y - X%*%beta)^8/(sigma2^6) - 24*nu^3*(y - X%*%beta)^10/(sigma2^7)))
  
  Hess[p + 2, p + 1,] <- Hess[p + 1, p + 2,]
  
  Hess[p + 2, p + 2,] <- ind*(-2*(10*(y - X%*%beta)^4/(sigma2^3)) + 
                                (12*(y - X%*%beta)^6/(sigma2^4) - 24*nu*(y - X%*%beta)^8/(sigma2^5) + 12*nu^2*(y - X%*%beta)^10/(sigma2^6)))
  
  return(Hess)
}

hess_H_score_Tukeys_regression_repar <- function(y, X, beta, sigma2, nu){
  
  n <- length(y)
  p <- length(beta)
  
  ind <- array(NA, dim = c(p + 2, p + 2, n))
  for(i in 1:(p + 2)){
    for(j in 1:(p + 2)){
      ind[i,j,] <- as.numeric((abs(y - X%*%beta) <= (sqrt(sigma2/nu))))
    }
  }
  return(ind * hess_H_score_Tukeys_regression_repar_noind(y, X, beta, sigma2, nu))
}

### Tukeys-loss Hessians - absapprox

#### Required functions

grad_msigmoid <- function(x, k){
  #return((-k*exp(k*x)) / ((1 + exp( k * x))^2))
  return(-exp(log(k) + k*x - 2*log1pexp( k * x)))
}

hessian_msigmoid <- function(x, k){
  return(exp(log(2) + 2*log(k) + 2*k*x - 3*log1pexp( k * x)) - exp(2*log(k) + k*x - 2*log1pexp( k * x)))
}

grad_abs_approx_regression <- function(y, X, beta, sigma2, k){
  
  p <- length(beta)
  n <- length(y)
  
  grad <- array(NA, dim = c(p + 1, n))
  for(i in 1:p){
    grad[i, ] <- -X[,i]*as.vector(y - X%*%beta)/(sigma2*sqrt(as.vector(y - X%*%beta)^2/sigma2 + 1/k))
  }
  grad[p + 1, ] <- -1/2*as.vector(y - X%*%beta)^2/(sigma2^2*sqrt(as.vector(y - X%*%beta)^2/sigma2 + 1/k))
  return(grad)
}

hessian_abs_approx_regression <- function(y, X, beta, sigma2, k){
  p <- length(beta)
  n <- length(y)
  Hess <- array(NA, dim = c(p + 1, p + 1, n))
  for(i in 1:p){
    for(j in 1:p){
      Hess[i, j, ] <- X[,i]*X[,j]*1/k/( sigma2*(as.vector(y - X%*%beta)^2/sigma2 + 1/k)^(3/2))
    }
    Hess[p + 1, i, ] <- X[,i]*as.vector(y - X%*%beta)*(as.vector(y - X%*%beta)^2/(2*sigma2) + 1/k)/(sigma2^2*(as.vector(y - X%*%beta)^2/sigma2 + 1/k)^(3/2))
    Hess[i, p + 1, ] <- Hess[p + 1, i, ]
  }
  Hess[p + 1, p + 1,] <- as.vector(y - X%*%beta)^2*(3*as.vector(y - X%*%beta)^2/(4*sigma2) + 1/k)/(sigma2^3*(as.vector(y - X%*%beta)^2/sigma2 + 1/k)^(3/2))
  return(Hess)
  
}


grad_approx_indicator_abs_less3_regression_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  
  n <- length(y)
  p <- length(beta)
  
  kappa <- 1/sqrt(nu)
  
  grad_abs_approx_eval <- grad_abs_approx_regression(y, matrix(X, nrow = n, ncol = p), beta, sigma2, k_abs)
  
  grad_sigmoid_eval <- - k_sigmoid * exp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa) - 2*log1pexp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa)))
  
  grad_msqrt_nu_m1 <-  1/(2*nu^(3/2))
  
  grad_beta <-  grad_abs_approx_eval[1:p,] * matrix(grad_sigmoid_eval, nrow = p, ncol = n, byrow = TRUE) 
  
  grad_sigma2 <- grad_abs_approx_eval[p + 1,] * matrix(grad_sigmoid_eval, nrow = 1, ncol = n, byrow = TRUE) 
  
  grad_nu <- grad_msqrt_nu_m1 * matrix(grad_sigmoid_eval, nrow = 1, ncol = n, byrow = TRUE) 
  
  grad <- array(NA, dim = c(p + 2, n))
  grad[1:p, ] <- grad_beta
  grad[p + 1, ] <- grad_sigma2
  grad[p + 2, ] <- grad_nu
  
  return(grad )
}

hessian_approx_indicator_abs_less3_regression_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  
  n <- length(y)
  p <- length(beta)
  
  kappa <- 1/sqrt(nu)
  
  grad_abs_approx_eval <- grad_abs_approx_regression(y, matrix(X, nrow = n, ncol = p), beta, sigma2, k_abs)
  hessian_abs_approx_eval <- hessian_abs_approx_regression(y, matrix(X, nrow = n, ncol = p), beta, sigma2, k_abs)
  
  grad_sigmoid_eval <- - k_sigmoid *exp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa) - 2*log1pexp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa)))
  
  hessian_sigmoid_eval <- exp(log(2) + 2*log(k_sigmoid) + 2*k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa) - 3*log1pexp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa))) - exp(2*log(k_sigmoid) + k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa) - 2*log1pexp(k_sigmoid*(abs_approx(as.vector(y - X%*%beta[1:p])/sqrt(sigma2), k_abs) - kappa)))
  
  grad_msqrt_nu_m1 <-  1/(2*nu^(3/2))
  hessian_msqrt_nu_m1 <-  -3/(4*nu^(5/2))
  
  Hess <- array(NA, dim = c(p + 2, p + 2, n))
  
  for(i in 1:(p + 1)){
    for(j in 1:(p + 1)){
      Hess[i, j, ] <- hessian_abs_approx_eval[i, j, ]*grad_sigmoid_eval + grad_abs_approx_eval[i, ]*grad_abs_approx_eval[j, ]*hessian_sigmoid_eval
    }
    Hess[p + 2, i, ] <- grad_msqrt_nu_m1*grad_abs_approx_eval[i, ]*hessian_sigmoid_eval
    Hess[i, p + 2, ] <- Hess[p + 2, i, ]
  }
  Hess[p + 2, p + 2, ] <- hessian_msqrt_nu_m1*grad_sigmoid_eval + grad_msqrt_nu_m1^2*hessian_sigmoid_eval
  
  
  return(Hess)
} 

hess_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  
  n <- length(y)
  p <- length(beta)
  
  H_score_Tukeys_regression_eval <- H_score_tukey_varThresh_noind(x = y, mu = X%*%beta, sigma2, 1/sqrt(nu))
  grad_H_score_Tukeys_regression_eval <- grad_H_score_Tukeys_regression_repar_noind(y, X, beta, sigma2, nu)
  hessian_H_score_Tukeys_regression_eval <- hess_H_score_Tukeys_regression_repar_noind(y, X, beta, sigma2, nu)
  
  tilde_ind <- approx_indicator_abs_less3((y - X%*%beta)/sqrt(sigma2), 1/sqrt(nu), k_abs, k_sigmoid)
  grad_tilde_ind <- grad_approx_indicator_abs_less3_regression_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid)
  hessian_tilde_ind <- hessian_approx_indicator_abs_less3_regression_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid)
  
  Hess <- array(NA, dim = c(p + 2, p + 2, n))
  for(i in 1:p){
    for(j in 1:p){
      Hess[i, j, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[i,j,] + H_score_Tukeys_regression_eval*hessian_tilde_ind[i,j,] + grad_tilde_ind[i,]*grad_H_score_Tukeys_regression_eval[j,] + grad_tilde_ind[j,]*grad_H_score_Tukeys_regression_eval[i,]
    }
    Hess[i, p + 1, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[i, p + 1, ] + H_score_Tukeys_regression_eval*hessian_tilde_ind[i, p + 1, ] + grad_tilde_ind[i,]*grad_H_score_Tukeys_regression_eval[p + 1,] + grad_tilde_ind[p + 1, ]*grad_H_score_Tukeys_regression_eval[i,]
    Hess[p + 1, i, ] <- Hess[i, p + 1, ]
    
    Hess[i, p + 2, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[i, p + 2, ] + H_score_Tukeys_regression_eval*hessian_tilde_ind[i, p + 2, ] + grad_tilde_ind[i,]*grad_H_score_Tukeys_regression_eval[p + 2,] + grad_tilde_ind[p + 2, ]*grad_H_score_Tukeys_regression_eval[i,]
    Hess[p + 2, i, ] <- Hess[i, p + 2, ]
  }
  
  Hess[p + 1, p + 1, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[p + 1, p + 1, ] + H_score_Tukeys_regression_eval*hessian_tilde_ind[p + 1, p + 1, ] + 2*grad_tilde_ind[p + 1, ]*grad_H_score_Tukeys_regression_eval[p + 1 ,]
  
  Hess[p + 1, p + 2, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[p + 1, p + 2, ] + H_score_Tukeys_regression_eval*hessian_tilde_ind[p + 1, p + 2, ] + grad_tilde_ind[p + 1,]*grad_H_score_Tukeys_regression_eval[p + 2,] + grad_tilde_ind[p + 2, ]*grad_H_score_Tukeys_regression_eval[p + 1,]
  Hess[p + 2, p + 1, ] <- Hess[p + 1, p + 2, ]
  
  Hess[p + 2, p + 2, ] <- tilde_ind*hessian_H_score_Tukeys_regression_eval[p + 2, p + 2, ] + H_score_Tukeys_regression_eval*hessian_tilde_ind[p + 2, p + 2, ] + 2*grad_tilde_ind[p + 2, ]*grad_H_score_Tukeys_regression_eval[p + 2 ,]
  
  return(Hess)
}

### Further functions for the Score Matching Information Criteria 

#### Gaussian model 

grad_squared_H_score_norm_regression <- function(y, X, beta, sigma2){
  
  n <- nrow(X)
  p <- length(beta)
  
  grad_beta <- matrix(NA, nrow = p, ncol = n)
  for(j in 1:p){
    grad_beta[j,] <- -2*X[,j]*(y - X%*%beta)/(sigma2^2)
  }
  grad_sigma2 <- 2/(sigma2^2) - 2*(y - X%*%beta)^2/(sigma2^3)
  
  grad_squared <- array(NA, dim = c(p + 1, p + 1, n))
  for(j in 1:p){
    for(k in 1:p){
      grad_squared[j, k,] <- grad_beta[j,]*grad_beta[k,]
    }
    grad_squared[j, p + 1, ] <- grad_beta[j,]*grad_sigma2
    grad_squared[p + 1, j, ] <- grad_squared[j, p + 1, ]
  }
  grad_squared[p + 1, p + 1, ] <- grad_sigma2^2
  return(grad_squared)
}


J_H_score_norm_regression <- function(y, X, beta, sigma2){
  apply(grad_squared_H_score_norm_regression(y, X, beta, sigma2), c(1, 2), mean)
}

K_H_score_norm_regression <- function(y, X, beta, sigma2){
  apply(hess_H_score_norm_regression(y, X, beta, sigma2), c(1, 2), mean)
}



#### Tukey's loss - regression - absapprox - repar

grad_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  
  n <- length(y)
  p <- length(beta)
  
  H_score_Tukeys_regression_eval <- H_score_tukey_varThresh_noind(x = y, mu = X%*%beta, sigma2, 1/sqrt(nu))
  grad_H_score_Tukeys_regression_eval <- grad_H_score_Tukeys_regression_repar_noind(y, X, beta, sigma2, nu)
  
  tilde_ind <- approx_indicator_abs_less3((y - X%*%beta)/sqrt(sigma2), 1/sqrt(nu), k_abs, k_sigmoid)
  grad_tilde_ind <- grad_approx_indicator_abs_less3_regression_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid)
  
  grad_beta <- array(NA, dim = c(p, n))
  for(i in 1:p){
    grad_beta[i,] <- H_score_Tukeys_regression_eval * grad_tilde_ind[i,] + tilde_ind * grad_H_score_Tukeys_regression_eval[i,]
  }
  
  grad_sigma2 <- H_score_Tukeys_regression_eval * grad_tilde_ind[p + 1,] + tilde_ind * grad_H_score_Tukeys_regression_eval[p + 1,]
  
  grad_nu <- H_score_Tukeys_regression_eval * grad_tilde_ind[p + 2,] + tilde_ind * grad_H_score_Tukeys_regression_eval[p + 2,]
  
  grad <- array(NA, dim = c(p + 2, n))
  grad[1:p,] <- grad_beta
  grad[p + 1,] <- grad_sigma2
  grad[p + 2,] <- grad_nu
  return(grad)
}

grad_squared_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  
  n <- length(y)
  p <- length(beta)
  
  grad_H_score_Tukeys_regression_absapprox_repar_eval <- grad_H_score_Tukeys_regression_absapprox_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid)
  
  grad_beta <- matrix(grad_H_score_Tukeys_regression_absapprox_repar_eval[1:p, ], nrow = p, ncol = n)
  
  grad_sigma2 <- grad_H_score_Tukeys_regression_absapprox_repar_eval[p + 1, ]
  
  grad_nu <- grad_H_score_Tukeys_regression_absapprox_repar_eval[p + 2, ]
  
  grad_squared <- array(NA, dim = c(p + 2, p + 2, n))
  for(i in 1:p){
    for(j in 1:p){
      grad_squared[i,j,] <- grad_beta[i,]*grad_beta[j,]
    }
    grad_squared[i, p + 1, ] <- grad_beta[i,]*grad_sigma2
    grad_squared[p + 1, i, ] <- grad_squared[i, p + 1, ]
    grad_squared[i, p + 2, ] <- grad_beta[i,]*grad_nu
    grad_squared[p + 2, i, ] <- grad_squared[i, p + 2, ]
  }
  grad_squared[p + 1, p + 1, ] <- grad_sigma2^2
  grad_squared[p + 1, p + 2, ] <- grad_sigma2*grad_nu
  grad_squared[p + 2, p + 1, ] <- grad_sigma2*grad_nu
  grad_squared[p + 2, p + 2, ] <- grad_nu^2
  return(grad_squared)
}

J_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  apply(grad_squared_H_score_Tukeys_regression_absapprox_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid), c(1, 2), mean)
}

K_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  apply(hess_H_score_Tukeys_regression_absapprox_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid), c(1, 2), mean)
}


SMIC_H_score_norm_regression <- function(y, X, beta, sigma2){
  insample_loss <- sum(H_score_norm(x = y, mu = X%*%beta, sigma2, w = 1))
  bias_adjustment <- sum(diag( J_H_score_norm_regression(y, X, beta, sigma2)%*%solve(K_H_score_norm_regression(y, X, beta, sigma2) )))
  
  return(insample_loss + bias_adjustment)
}


SMIC_H_score_Tukeys_regression_absapprox_repar <- function(y, X, beta, sigma2, nu, k_abs, k_sigmoid){
  insample_loss <- sum(H_score_tukey_varThresh_absapprox(x = y, mu = X%*%beta, sigma2, c = 1/sqrt(nu), k_abs, k_sigmoid))
  bias_adjustment <- sum(diag(J_H_score_Tukeys_regression_absapprox_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid)%*%solve(K_H_score_Tukeys_regression_absapprox_repar(y, X, beta, sigma2, nu, k_abs, k_sigmoid) )))
  
  return(insample_loss + bias_adjustment)
}

