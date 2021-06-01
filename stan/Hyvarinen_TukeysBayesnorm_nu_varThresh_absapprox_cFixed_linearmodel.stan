
functions {
   
   // Here nu = 1/kappa^2
   
   real sigmoid(real x, real k){
      // return (1/(1+exp( - k * x)));
      return exp( - log1p_exp( - k * x));
   }
   
   real abs_approx(real x, real k){
      return sqrt(x^2 + 1/k);
   }
   
   real approx_indicator_abs_less3(real x, real c, real k){
      return sigmoid( - (abs_approx(x, k) - c), k);
   }
   
   real bdown_penalty_lpdf(real bdown, real penalty){
      return penalty * sigmoid(bdown, -100);
   }
   
   real Tukeys_score_norm (real x, real mu, real sigma2, real nu, real w){
     int ind;
     ind = (((x - mu)^2/sigma2) <= (1/nu));
     
     return w*(ind*(((x-mu)^2/(2*sigma2) - (nu*(x-mu)^4)/(2*sigma2^2) + (nu^2*(x-mu)^6)/(6*sigma2^3))) + (1-ind)*(1/(6*nu)));
   }
   
   
   real grad_Tukeys_score_norm (real x, real mu, real sigma2, real nu, real w){
     real ind;
     ind = approx_indicator_abs_less3((x - mu)/sqrt(sigma2), 1/sqrt(nu), 100);
     
     return ind*(-w*((x-mu)/sigma2 - 2*(nu*(x-mu)^3)/(sigma2^2) + (nu^2*(x-mu)^5)/(sigma2^3)));
   }

   real Laplacian_Tukeys_score_norm(real x, real mu, real sigma2, real nu, real w){
     real ind;
     ind = approx_indicator_abs_less3((x - mu)/sqrt(sigma2), 1/sqrt(nu), 100);
     
     return ind*(-w*(1/sigma2 - 6*(nu*(x-mu)^2)/(sigma2^2) + 5*(nu^2*(x-mu)^4)/(sigma2^3)));
   }

   real H_score_Tukeys_lpdf (real x, real mu, real sigma2, real nu, real w){
      
      return 2*Laplacian_Tukeys_score_norm(x, mu, sigma2, nu, w) + (grad_Tukeys_score_norm(x, mu, sigma2, nu, w))^2;
   }


}


data {
   
   int<lower=0> n;
   int<lower=0> p;
   matrix[n,1] y;
   matrix[n,p] X;
   real mu_beta;
   real<lower=0> beta_s;
   real<lower=0> sig_p1;
   real<lower=0> sig_p2;
   real<lower=0> w;
   real<lower=0> nu;
}

parameters 
{
   
   vector[p] beta;
   real<lower = 0> sigma2;

}

transformed parameters
{
   matrix[n,1] lin_pred;
   lin_pred[,1] = X*beta;

}


model {
   

   real tukey_loss_thresh;
   real bdown_threshold;
   tukey_loss_thresh = 0;
   for(i in 1:n){
      tukey_loss_thresh  += Tukeys_score_norm (y[i,1], lin_pred[i,1], sigma2, nu, 1);
   }
   
   
   bdown_threshold = (1/(n/2.0 - p + 1))*(((n/2.0 - p + 1)) - (tukey_loss_thresh/(1/(6.0*nu))));
   
   target += inv_gamma_lpdf(sigma2 | sig_p1, sig_p2);
   target += normal_lpdf(beta | mu_beta, sqrt(sigma2*beta_s));
   
   target += bdown_penalty_lpdf(bdown_threshold | - 10000);

   for(i in 1:n){
     target += -w*H_score_Tukeys_lpdf(y[i,1]| lin_pred[i,1], sigma2, nu, 1);
   }
}

