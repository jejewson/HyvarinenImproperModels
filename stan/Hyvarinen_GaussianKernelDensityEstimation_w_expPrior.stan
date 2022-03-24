
functions {
   
   matrix Gaussian_K_matrix(vector data_total){
      int n = num_elements(data_total) ;
      matrix[n, n] K;
      for(j in 1:n){
         K[,j] = (data_total - data_total[j]);
      }
      return K;
   }
   
   real H_score_KDE_Gaussian_simplified_lpdf (real sigma2, row_vector K_i, real w){
      int n;
      real numerator1;
      real denominator1;
      real numerator2;
      n = num_elements(K_i);
      numerator1 =  sum(square(K_i)/(n*sqrt(2*pi())*sigma2^(2.5)) .* exp( - square(K_i)/(2.0*sigma2))) - 
                         sum(1/(n*sqrt(2*pi())*sigma2^(1.5))*exp( - square(K_i)/(2.0*sigma2)));
      denominator1 = sum(1/(n*sqrt(2*pi())*sigma2^(0.5))*exp( - square(K_i)/(2.0*sigma2)));
      numerator2 =  sum((K_i)/(n*sqrt(2*pi())*sigma2^(1.5)) .* exp( - square(K_i)/(2.0*sigma2)));
      numerator1 -= - 1/(n*sqrt(2*pi())*sigma2^(1.5));// for the extra numerator 1 when data_i = data_total[i]
      return 2.0*w*numerator1/denominator1 - (2.0*w - w^2)*numerator2^2/denominator1^2;
   }
   

}


data {
   
   int<lower=0> n;
   matrix[n,1] y;
   real<lower=0> sig_p1;
   real<lower=0> sig_p2;
   real<lower=0> omega_p1;
   real<lower=0> w;
   real<lower=0> sigma2_lower;
}

parameters 
{
   
   real<lower = sigma2_lower> sigma2;
   real<lower = 0> omega;
}

model {
   
   matrix[n, n] K;
   K = Gaussian_K_matrix(y[,1]);
   
   
   target += inv_gamma_lpdf(sigma2 | sig_p1, sig_p2);
   target += exponential_lpdf(omega | omega_p1);

   for(i in 1:n){
     target += - w * H_score_KDE_Gaussian_simplified_lpdf(sigma2 | K[i, ], omega);
   }
}

