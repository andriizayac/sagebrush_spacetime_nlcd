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
  real<lower=0,upper=1> r;          // global growt rate
  real<lower=0,upper=1> r_mu;       // avg growth 
  real<lower=0> r_sigma;    // growth variance
  vector<lower=0>[N] k;     // pixel-level carrying capacity
  real<lower=0> k_mu;       // average (global) k
  real<lower=0> k_sigma;    // k variance 
  real k_phi;               // change in post-burn caryying capacity
  vector<lower=0>[T] eta;   // observation error
  vector[T-1] gamma;        // observation error
  
  vector<lower=0,upper=1>[N] z_init;          // initial state vector (non-dim in k)
  vector<lower=0,upper=1>[N_pred] z_pred;     // initial state vector (non-dim in k)
  //vector<lower=6>[N_pred] k;
} transformed parameters {
  vector<lower=0>[N] Z[T]; // latent states
  Z[1] = z_init; 
  for(t in 2:T)
      Z[t] = Z[t-1] + r*Z[t-1] - r*(Z[t-1] .* Z[t-1]);//+ gamma[t-1]; // )* beta[t]
} model {
  //priors
  alpha ~ normal(0,1);
  r_mu ~ normal(.5,1)T[0,1];
  r_sigma ~ exponential(5);
  r ~ normal(r_mu,r_sigma);
  k_mu ~ normal(0, 5);
  k_sigma ~ exponential(1);
  k_phi ~ normal(0,.1)T[-1,1];
  eta ~ exponential(1); 
  gamma ~ normal(0,.01);
  for(i in 1:N) z_init[i] ~ normal(0,.5)T[0,1];
  for(i in 1:N_pred) z_pred[i] ~ normal(0,.5)T[0,1];

  
  k ~ normal(k_mu + k_phi*y_k_prior_sc,k_sigma);

  y0_tr ~ normal(z_init,eta[1]);
  y0_pred ~ normal(z_pred,.1);
   
   for(t in 1:T){
      y_mat[t] ~ normal(alpha + Z[t] .* k,eta[t]);
   }
} generated quantities {
  matrix[T_pred,N_pred] y_pred01;
  matrix[T_pred,N_pred] y_pred;
  //vector<lower=0>[N_pred] k_pred = mean(y_k_pred) + k_phi*y_k_pred_sc;

  for(i in 1:N_pred)  y_pred01[1,i] = z_pred[i];
  
  for(t in 2:T_pred){
    for(i in 1:N_pred)
      y_pred01[t,i] = (y_pred01[t-1,i]+ r*y_pred01[t-1,i] - r*y_pred01[t-1,i]^2);
    }
  for(t in 1:T_pred) 
    y_pred[t] = alpha+y_pred01[t] .* (mean(y_k_pred) + k_phi*y_k_pred_sc');
    
}
