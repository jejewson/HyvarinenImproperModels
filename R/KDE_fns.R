
### Laplace Approximation Functions

log_laplace_approximation_marg_lik <- function(log_P_star, A, p){
  return(log_P_star + p/2*log(2*pi) - 1/2*log(det(A)))
}

### Hyvarinen score functions

f_KDE_Gaussian <- function(x, data, sigma2, w){
  n <- length(data)
  return((1/n*sum(dnorm(x, mean = data, sd = sqrt(sigma2))))^w)
}

f_KDE_Gaussian_vect <- function(x, data, sigma2, w){
  m <- length(x)
  n <- length(data)
  out <- rep(NA, m)
  for(j in 1:m){
    out[j] <- (1/n*sum(dnorm(x[j], mean = data, sd = sqrt(sigma2))))^w
  }
  
  return(out)
}

log_f_KDE_Gaussian <- function(x, data, sigma2, w){
  n <- length(data)
  return(w*(log(1/n) + logSumExp(dnorm(x, mean = data, sd = sqrt(sigma2), log = TRUE))))
}

## In-sample training

f_KDE_Gaussian_insamp <- function(i, data, sigma2, w){
  n <- length(data)
  return((1/n*sum(dnorm(data[i], mean = data, sd = sqrt(sigma2))))^w)
}

log_f_KDE_Gaussian_insamp <- function(data_i, data_m, sigma2, w){
  n <- length(data_m) + 1
  return(w*(log(1/n) + logSumExp(dnorm(data_i, mean = c(data_i, data_m), sd = sqrt(sigma2), log = TRUE))))
}

grad_log_f_KDE_Gaussian_insamp <- function(data_i, data_m, sigma2, w){
  n <- length(data_m) + 1
  return(w*(-sum((data_i - c(data_i, data_m))/(n*sqrt(2*pi)*sigma2^(3/2))*exp(-(data_i - c(data_i, data_m))^2/(2*sigma2)))/(sum(1/(n*sqrt(2*pi)*sigma2^(1/2))*exp(-(data_i - c(data_i, data_m))^2/(2*sigma2))))))
}

Laplacian_log_f_KDE_Gaussian_insamp <- function(data_i, data_m, sigma2, w){
  n <- length(data_m) + 1
  return(w*(
    (sum((data_i - data_m)^2/(n*sqrt(2*pi)*sigma2^(5/2))*exp(-(data_i - data_m)^2/(2*sigma2))) - sum(1/(n*sqrt(2*pi)*sigma2^(3/2))*exp(-(data_i - data_m)^2/(2*sigma2))))/(1/n*sum(dnorm(data_i, mean = c(data_i, data_m), sd = sqrt(sigma2))))
    -(sum((data_i - data_m)/(n*sqrt(2*pi)*sigma2^(3/2))*exp(-(data_i - data_m)^2/(2*sigma2))))^2/((1/n*sum(dnorm(data_i, mean = c(data_i, data_m), sd = sqrt(sigma2))))^2)
  ))
}

H_KDE_Gaussian_insamp <- function(data_i, data_m, sigma2, w){
  return(2*Laplacian_log_f_KDE_Gaussian_insamp(data_i, data_m, sigma2, w) + (grad_log_f_KDE_Gaussian_insamp(data_i, data_m, sigma2, w))^2)
}

### Hyvarinen-score sample functions

sample_log_f_KDE_Gaussian_insamp <- function(data, sigma2, w){
  n <- length(data)
  out <- 0
  for(i in 1:n){
    out <- out + log_f_KDE_Gaussian_insamp(data[i], data[-i], sigma2, w)
  }
  return(out)
}

sample_H_KDE_Gaussian_insamp <- function(data, sigma2, w){
  n <- length(data)
  out <- 0
  for(i in 1:n){
    out <- out + H_KDE_Gaussian_insamp(data[i], data[-i], sigma2, w)
  }
  return(out)
}

### Dirichlet Process Mixture Model

dnormDPmixture_MAP <- function(x, weights, means, sds){
  
  num_components <- length(weights)
  num_obs <- length(x)
  dens <- rep(0, num_obs)
  for(j in 1:num_components){
    dens <- dens + weights[j]*dnorm(x, mean = means[j], sd = sds[j])
  }
  return(dens)
  
}

### Gaussian Mixture Model Density

dnorm_mixture <- function(x, mus, sigma2s, ws){
  num_comp <- length(mus)
  num_obs <- length(x)
  dens <- rep(0, num_obs)
  for(j in 1:num_comp){
    dens <- dens + ws[j]*dnorm(x, mus[j], sqrt(sigma2s[j]))
  }
  return(dens)
}

### Numerically estimating Fisher's Divergence functions

FisherDivergence <- function(g, grad_log_g, grad_log_f, lower = -Inf, upper = Inf){
  
  FisherDivergence_integrand <- function(x){((grad_log_g(x) - grad_log_f(x))^2)*g(x)}
  FisherDivergence_integral <- integrate(f = FisherDivergence_integrand, lower, upper, subdivisions = 1000, rel.tol = .Machine$double.eps^0.5)# defaults , subdivisions = 100, rel.tol = .Machine$double.eps^0.25)
  
  return(FisherDivergence_integral)
}


grad_dnorm <- function(x, mu, sigma2){
  return(-(x-mu)/(sqrt(2*pi)*sigma2^(3/2))*exp(-(x-mu)^2/(2*sigma2)))
}

grad_log_dnorm_mixture <- function(x, mus, sigma2s, ws){
  num_comp <- length(mus)
  num_obs <- length(x)
  dens <- rep(0, num_obs)
  grad <- rep(0, num_obs)
  for(j in 1:num_comp){
    dens <- dens + ws[j]*dnorm(x, mus[j], sqrt(sigma2s[j]))
    grad <- grad + ws[j]*grad_dnorm(x, mu = mus[j], sigma2 = sigma2s[j])
  }
  return(grad/dens)
}

grad_dnorm_std <- function(x, xbar, s, mu, sigma2){
  ## grad dnorm s*x + xbar
  return(-s^2*(s*x + xbar - mu)/(sqrt(2*pi)*sigma2^(3/2))*exp(-(s*x + xbar - mu)^2/(2*sigma2)))
}

grad_log_dnorm_mixture_std <- function(x, xbar, s, mus, sigma2s, ws){
  num_comp <- length(mus)
  num_obs <- length(x)
  dens <- rep(0, num_obs)
  grad <- rep(0, num_obs)
  for(j in 1:num_comp){
    dens <- dens + ws[j]*s*dnorm(s*x + xbar, mus[j], sqrt(sigma2s[j]))
    grad <- grad + ws[j]*grad_dnorm_std(x, xbar, s, mu = mus[j], sigma2 = sigma2s[j])
  }
  return(grad/dens)
}


grad_log_f_KDE_Gaussian_vect <- function(x, data, sigma2, w){
  m <- length(x)
  n <- length(data)
  out <- rep(NA, m)
  for(j in 1:m){
    out[j] <- w*(-sum((x[j] - data)/(n*sqrt(2*pi)*sigma2^(3/2))*exp(-(x[j] - data)^2/(2*sigma2)))/(sum(1/(n*sqrt(2*pi)*sigma2^(1/2))*exp(-(x[j] - data)^2/(2*sigma2)))))
  }  
  return(out)
}

