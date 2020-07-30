data{
  int N;
  int T;
  matrix[N,T] y_mat;
} 
transformed data{
  matrix[T,N] y = y_mat';
}
parameters{
  vector[T] beta;
  real<lower=0> eta;
}
model{
  beta ~ normal(0,10);
  eta ~ exponential(1);
  for(t in 1:T)
    y[t] ~ normal(beta[t],eta);
}
