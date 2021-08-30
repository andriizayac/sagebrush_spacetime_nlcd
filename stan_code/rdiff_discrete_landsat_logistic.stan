data {
  int N;  // number of observations in the matrix
  int n;
  int T;  // number of time steps in a series
  vector[N] y_mat[T]; // response, cover ts
  vector[N] y0;
  vector[n] y0pred;
  vector[N] y_k_prior;
  vector[n] y_k_pred;
} transformed data {
  vector[N] y_k_prior_sc = (y_k_prior - mean(y_k_prior)) ; //./ 2*sd(y_k_prior)
  vector[n] y_k_pred_sc = (y_k_pred - mean(y_k_pred)) ; //./ 2*sd(y_k_prior)
} parameters {
  vector<lower=0>[N] r;          // pixel-level growt rate
  real<lower=0,upper=1> r_mu;       // avg growth 
  real<lower=0> r_sigma;    // growth variance
  vector<lower=0>[N] k;     // pixel-level carrying capacity
  real<lower=0> k_mu;       // average (global) k
  real<lower=0> k_sigma;    // k variance 
  real<lower=0,upper=1> k_phi;               // change in post-burn caryying capacity
  vector<lower=0>[T-1] eta;   // observation error
  
} transformed parameters {
} model {
  // priors: demographic parameters
  r_mu ~ normal(.5, .5);
  r_sigma ~ exponential(5);
  r ~ normal(r_mu, r_sigma);
  k_mu ~ normal(0, 5);
  k_sigma ~ exponential(1);
  k_phi ~ normal(0, .1);
  eta ~ exponential(1); 

  k ~ normal(k_mu + k_phi * y_k_prior_sc, k_sigma);

   for(t in 2:T){
     for(i in 1:N){
      y_mat[i, t] ~ normal(y_mat[i,t-1] + r[i] * y_mat[i,t-1] - r[i] * (y_mat[i,t-1]^2) / k[i], eta[t]);
   }
  }
} generated quantities {
  /*
  matrix[T,n] y_pred01;
  matrix[T,n] y_pred;

  for(i in 1:n)  y_pred01[1,i] = z0[i];
  
  for(t in 2:T)
    for(i in 1:n)
      y_pred01[t,i] = (y_pred01[t-1,i]+ r_mu * y_pred01[t-1,i] - r_mu * y_pred01[t-1,i]^2);
    
  for(t in 1:T) 
    y_pred[t] = alpha + y_pred01[t] .* (k_mu  + k_phi * y_k_pred_sc');
  */  
}
