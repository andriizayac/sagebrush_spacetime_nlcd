functions{
  real harm_mean(vector D, int N){
    real harmonic_mean;
    vector[2] n = to_vector({N,1});
    harmonic_mean = n[1]/sum(1 ./ D);
    return harmonic_mean;
  }
} 
data{
  int N;  // number of observations in the matrix
  int T;  // number of time steps in a series
  vector[N] y_mat[T]; // response, cover ts
  vector[N] y0_tr;
  vector[N] y_k_prior;
  int NN[N, 4];
} 
transformed data{
  real dx = 30;
  real dy = 30;
  real dt = 1;
  real dhat = dt/(dx^2 + dy^2);
  real dhatID =  -4 * dt/(dx^2 + dy^2); 
  // === spatial prior
  // vector[N] zeros;
  // matrix<lower = 0>[N, N] D;
  // === spatial prior
  // {
  //   vector[N] W_rowsums;
  //   for (i in 1:N) {
  //     W_rowsums[i] = sum(W[i, ]);
  //   }
  //   D = diag_matrix(W_rowsums);
  // }
  // zeros = rep_vector(0, N);
} parameters {
  //real alpha;
  real<lower=0> r;          // global growt rate
  real<lower=0> r_mu;       // avg growth 
  real<lower=0> r_sigma;    // growth variance
  real<lower=0> k_mu;       // average (global) k
  real<lower=0> k_sigma;    // k variance 
  real k_phi;               // change in post-burn caryying capacity
  real<lower=0> gamma;        // process error
  
  vector<lower=0>[N] k;     // pixel-level carrying capacity
  vector<lower=0>[N] Delta;
  //real<lower=0> delta;
  vector<lower=0>[T] eta;   // observation error

  vector<lower=0,upper=1>[N] z_init;          // initial state vector (non-dim in k)
  vector<lower=0>[N] Z[T+1]; // latent states
  // === spatial prior
  // vector[n] phi;
  // real<lower = 0> tau;
  // real<lower = 0, upper = 1> alpha;
} model {
  // === priors: intercepts
  //alpha ~ normal(0, 1);
  // === priors: demographic parameters
  r_mu ~ normal(.5, 1);
  r_sigma ~ exponential(5);
  r ~ normal(r_mu, r_sigma);
  k_mu ~ normal(0, 5);
  k_sigma ~ exponential(1);
  k_phi ~ normal(0, 1)T[-1, 1];
  // === variances
  eta ~ exponential(1); 
  gamma ~ exponential(5);
  // === spatial prior
  Delta ~ gamma(3, .1);
  // delat ~ normal(50, 5);
  // phi ~ multi_normal_prec(zeros, tau * (D - alpha * W));
  // Delta ~ normal(delta + phi, 5);
  
  for(i in 1:N) z_init[i] ~ lognormal(0, 1);
  y0_tr ~ normal(z_init, gamma);
  Z[1] ~ normal(z_init, gamma); 
  
  for(t in 2:T+1){
    for(j in 1:N){
      Z[t, j] ~ normal(Z[t-1, j] +
        Z[t-1, j] * Delta[j] * dhatID +
        Z[t-1, NN[j, 1]] * Delta[j] * dhat +
        Z[t-1, NN[j, 2]] * Delta[j] * dhat +
        Z[t-1, NN[j, 3]] * Delta[j] * dhat +
        Z[t-1, NN[j, 4]] * Delta[j] * dhat +
        r * Z[t-1, j] - r * (Z[t-1, j] * Z[t-1, j]), gamma);
    }
  }
  k ~ normal(k_mu + k_phi * y_k_prior, k_sigma);
  
   for(t in 2:T+1){
      y_mat[t-1] ~ normal(Z[t] .* k, eta[t-1]);
   }
}
