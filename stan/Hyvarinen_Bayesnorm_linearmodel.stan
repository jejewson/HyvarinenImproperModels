
functions {
   
   real grad_score_norm (real y, row_vector X, vector beta, real sigma2, real w){

     return (-w*(y-X*beta)/sigma2);
   }

   real Laplacian_score_norm(real y, row_vector X, vector beta, real sigma2, real w){
     
     return -w/sigma2;
   }

   real H_score_norm_lpdf (real y, row_vector X, vector beta, real sigma2, real w ){
      
      return 2*Laplacian_score_norm(y, X, beta, sigma2, w) + (grad_score_norm(y, X, beta, sigma2, w))^2;
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
}

parameters 
{
   
   vector[p] beta;
   real<lower=0> sigma2;

}

transformed parameters
{
   matrix[n,1] lin_pred;
   lin_pred[,1] = X*beta;
}


model {

   target += inv_gamma_lpdf(sigma2 | sig_p1, sig_p2);
   target += normal_lpdf(beta | mu_beta,sqrt(sigma2*beta_s));

   for(i in 1:n){
     target += -w*H_score_norm_lpdf(y[i,1]| X[i,], beta, sigma2, 1);
   }
}

