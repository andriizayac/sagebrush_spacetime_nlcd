data {
  int N;  // number of observations in the matrix
  int T;  // number of time steps in a series
  int N_pred;
  int T_pred;
  vector[N] y_mat[T]; // response, cover ts
  vector[N] y0_tr;
  vector[N_pred] y0_pred;
  vector[N] y_k_prior;
  vector[N_pred] y_k_pred;
} transformed data {
  vector[N] y_k_prior_sc = (y_k_prior - mean(y_k_prior)) ; //./ 2*sd(y_k_prior)
  vector[N_pred] y_k_pred_sc = (y_k_pred - mean(y_k_pred)) ; //./ 2*sd(y_k_pred)
} parameters {
  real alpha;
  real alpha_pred;
  real<lower=0> sigma_alpha_pred;
  real<lower=0,upper=1> r;          // global growt rate
  real<lower=0,upper=1> r_mu;       // avg growth 
  real<lower=0> r_sigma;    // growth variance
  vector<lower=0>[N] k;     // pixel-level carrying capacity
  real<lower=0> k_mu;       // average (global) k
  real<lower=0> k_sigma;    // k variance 
  real k_phi;               // change in post-burn caryying capacity
  vector<lower=0>[T] eta;   // observation error
  real<lower=0> gamma;        // observation error
  
  vector<lower=0,upper=1>[N] z_init;          // initial state vector (non-dim in k)
  vector<lower=0,upper=1>[N_pred] z_pred;     // initial state vector (non-dim in k)
  
  vector<lower=0>[N] Z[T+1]; // latent states
} transformed parameters {

} model {
  // priors: intercempts
  alpha ~ normal(0, 1);
  alpha_pred ~ normal(0, 10);
  sigma_alpha_pred ~ exponential(1);
  // priors: demographic parameters
  r_mu ~ normal(.5, 1)T[0, 1];
  r_sigma ~ exponential(5);
  r ~ normal(r_mu, r_sigma);
  k_mu ~ normal(0, 5);
  k_sigma ~ exponential(1);
  k_phi ~ normal(0, .1)T[-1, 1];
  eta ~ exponential(1); 
  gamma ~ exponential(1);
  
  for(i in 1:N) z_init[i] ~ lognormal(0, 1);
  for(i in 1:N_pred) z_pred[i] ~ lognormal(0, 1);
  
  Z[1] ~ normal(z_init,gamma); 
  for(t in 2:T+1)
      Z[t] ~ normal(Z[t-1] + r*Z[t-1] - r*(Z[t-1] .* Z[t-1]), gamma);
  
  k ~ normal(k_mu + k_phi * y_k_prior_sc, k_sigma);
  y_k_pred ~ normal(alpha_pred, sigma_alpha_pred);

  y0_tr ~ normal(z_init, eta[1]);
  y0_pred ~ normal(z_pred, .1);
   
   for(t in 2:T+1){
      y_mat[t-1] ~ normal(alpha + Z[t] .* k, eta[t]);
   }
} generated quantities {
  matrix[T_pred, N_pred] y_pred01;
  matrix[T_pred, N_pred] y_pred;
  //vector<lower=0>[N_pred] k_pred = mean(y_k_pred) + k_phi*y_k_pred_sc;

  for(i in 1:N_pred)  y_pred01[1, i] = z_pred[i];
  
  for(t in 2:T_pred){
    for(i in 1:N_pred)
      y_pred01[t, i] = (y_pred01[t-1, i]+ r * y_pred01[t-1, i] - r * y_pred01[t-1, i]^2);
    }
  for(t in 1:T_pred) 
    y_pred[t] = alpha + y_pred01[t] .* (alpha_pred + k_phi * y_k_pred_sc');
    
}
