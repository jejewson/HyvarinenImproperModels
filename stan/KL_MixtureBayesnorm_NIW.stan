// Mixture of gaussians using the local NormalInverseWishart prior 
// used in mombf's 'bfnormmix' function

// Note the univariate case of IW(Sigma_j; nu0, S0) is inverse-gamma distribution 
// with alpha0 = nu0/2 and \beta0 = S0/2 

functions {
   
     // Creates one simplex
  vector simplex_create(vector theta_raw, int m){
    vector[m+1] theta_simplex;
    real stick_len = 1;
    real prop;
    for(j in 1:m){
      prop = inv_logit(theta_raw[j] - log(m - j + 1));
      theta_simplex[j] = stick_len * prop;
      stick_len = stick_len - theta_simplex[j];
    }
    theta_simplex[m + 1] = stick_len;
    
    return theta_simplex;
  }
  
  // To correctly specify we also need the log absolute Jacobian determinant of the simplex_create
  real simplex_create_lj(vector theta_raw, int m){
    real lj = 0;
    real stick_len = 1;
    real adj_theta_raw_j;
    for (j in 1:m) {
      adj_theta_raw_j = theta_raw[j] - log(m - j + 1);
      lj = lj + log(stick_len) - log1p_exp(-adj_theta_raw_j) - log1p_exp(adj_theta_raw_j);
      stick_len = stick_len * (1.0 - inv_logit(adj_theta_raw_j));
    }
    
    return lj;
  }
   
   real norm_mix_lpdf (real y, int K, vector mu, vector sigma2, vector omega){
      real log_lik = negative_infinity();
      for(k in 1:K){
         log_lik = log_sum_exp(log_lik, log(omega[k]) + normal_lpdf(y | mu[k], sqrt(sigma2[k])));
      }
     return log_lik;
   }


}


data {
   
   int<lower=0> n;
   int<lower=1> K;
   matrix[n,1] y;
   vector[K] mu_0;
   real<lower=0> kappa; // g in mombf
   real<lower=0> nu_0;
   real<lower=0> S_0;
   real<lower=0> alpha_0; // q.niw in mombf

}

parameters 
{
   
   vector[K] mu;
   vector<lower=0>[K] sigma2;
   vector[K - 1] omega_raw; // Each K simplex only has K - 1 degrees of freedom

}

model {
   
   vector[K] omega_simplex;
   // First we turn our raw omega's into simplexes
   omega_simplex = simplex_create(omega_raw, K - 1);


   target += dirichlet_lpdf(omega_simplex | rep_vector(alpha_0, K)) + 
               simplex_create_lj(omega_raw, K - 1);
   
   target += inv_gamma_lpdf(sigma2 | 0.5*nu_0, 0.5*S_0);
   target += normal_lpdf(mu | mu_0, sqrt(kappa*sigma2));
   
   for(i in 1:n){
     target += norm_mix_lpdf(y[i, 1] | K, mu , sigma2, omega_simplex);
   }
}

generated quantities {
   vector[K] omega_simplex;
   omega_simplex = simplex_create(omega_raw, K - 1);
   
}

