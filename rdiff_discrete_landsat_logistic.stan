data {
  int N;  // number of observations in the matrix
  int T;  // number of time steps in a series
  vector[N] y_mat[T]; // response, cover ts
  vector[N] y0_tr;
} transformed data {
} parameters {
  vector<lower=0> r;          // pixel-level growt rate
  //real<lower=0,upper=1> r_mu;       // avg growth 
  //real<lower=0> r_sigma;    // growth variance
  vector<lower=0> k;     // pixel-level carrying capacity
  //real<lower=0> k_mu;       // average (global) k
  //real<lower=0> k_sigma;    // k variance 
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
}
