functions {
  real harm_mean(vector D, int N){
    real harmonic_mean;
    vector[2] n = to_vector({N,1});
    harmonic_mean = n[1]/sum(1 ./ D);
    return harmonic_mean;
  }
  matrix H_mat(vector Dvec, real r, real dx, real dy, int[,] NN, real[,] NN01, int N){
    matrix[N,N] H;
    H = rep_matrix(0.0,N,N);
    for(i in 1:N){
      H[i,i] = Dvec[i]*(1/dx^2 + 1/dy^2);
      if(sum(NN01[i]) == 4){
        if(NN[i,1] > 0){H[i,NN[i,1]] = (Dvec[i] - (Dvec[NN[i,2]] - Dvec[NN[i,1]])/4)/(2*dx^2);}
        if(NN[i,2] > 0){H[i,NN[i,2]] = (Dvec[i] + (Dvec[NN[i,2]] - Dvec[NN[i,1]])/4)/(2*dx^2);}
        if(NN[i,3] > 0){H[i,NN[i,3]] = (Dvec[i] + (Dvec[NN[i,4]] - Dvec[NN[i,3]])/4)/(2*dy^2);}
        if(NN[i,4] > 0){H[i,NN[i,4]] = (Dvec[i] - (Dvec[NN[i,4]] - Dvec[NN[i,3]])/4)/(2*dy^2);}
    }}
      return H;
  }
  matrix G_mat(vector u, real dx, real dy,  int[,] NNind, int[,] NN, int N){
    matrix[N,5] G;
    matrix[N,4] ind;
    G = rep_matrix(0.0,N,5);
    for(i in 1:N){
      for(j in 1:4){
        ind[i,j] = NN[i,j] > 0 ? 1 : 0;
      }
    }
    for(i in 1:N){
      if(sum(ind[i])==4){
        G[i,3] = (u[NN[i,2]] - 2*u[i] + u[NN[i,1]])/(2*pow(dx,2)) + (u[NN[i,4]] - 2*u[i] + u[NN[i,3]])/(2*pow(dy,2));
        if(NN[i,1]>0){G[i,NNind[i,1]] = -(u[NN[i,2]] - u[NN[i,1]])/(4*pow(dx,2));}
        if(NN[i,2]>0){G[i,NNind[i,2]] = (u[NN[i,2]] - u[NN[i,1]])/(4*pow(dx,2));}
        if(NN[i,3]>0){G[i,NNind[i,3]+1] = -(u[NN[i,4]] - u[NN[i,3]])/(4*pow(dy,2));}
        if(NN[i,4]>0){G[i,NNind[i,4]+1] = (u[NN[i,4]] - u[NN[i,3]])/(4*pow(dy,2));}
    }}
      return G;
  }
  row_vector H_vec(vector D, real r, real dx, real dy, real[] NN01){
    row_vector[5] h;
    // h[3] = 1 + r - D[1]*(1/dx^2+1/dy^2); // central patch
    // h[1] = (-D[3] + 2*D[1] + D[2])/(4*dx^2); // left neighbor
    // h[2] = (D[3] + 2*D[1] - D[2])/(4*dx^2); // right neighbor
    // h[4] = (-D[5] + 2*D[1] + D[4])/(4*dy^2); // bottom neighbor
    // h[5] = (D[5] + 2*D[1] - D[4])/(4*dy^2); // top naighbor
    h[3] = 1-2* D[1]*(1/dx^2+1/dy^2); // central patch
    h[1] = (D[1] - (D[3] - D[2])/4)/dx^2; // left neighbor
    h[2] = (D[1] + (D[3] - D[2])/4)/dx^2; // right neighbor
    h[4] = (D[1] + (D[5] - D[4])/4)/dy^2; // bottom neighbor
    h[5] = (D[1] - (D[5] - D[4])/4)/dy^2; // top naighbor
    //for (i in 1:5)
      //h[i] = NN01[i] > 0 ? h[i] : 0.0;
    return h;
  }
  row_vector G_vec(vector u, real dx, real dy, real[] NN01){
    row_vector[5] g;
    g[1] = (u[3] - 2*u[1] + u[2])/(2*dx^2) + (u[5] - 2*u[1] + u[4])/(2*dy^2);
    g[2] = -(u[3] - u[2])/(4*dx^2);
    g[3] = (u[3] - u[2])/(4*dx^2);
    g[4] = -(u[5] - u[4])/(4*dy^2);
    g[5] = (u[5] - u[4])/(4*dy^2);
    for (i in 1:5)
      g[i] = NN01[i] > 0 ? g[i] : 0.0;
    return g;
  }
  vector delta_vec(vector D, int[] NNind, real[] NN01){
    vector[5] d;
    d = D[NNind];
    for (i in 1:5)
      d[i] = NN01[i] > 0 ? d[i] : 0.0;
    return d;
  }
} data {
  int N;  // number of observations in the matrix
  int T;  // number of time steps in a series
  vector[N] y_mat[T]; // response, cover ts
  int NN[N,4];
  int NNind[N,4];
  real NN01[N,4];
  int N_tr;  // number of observations in the matrix
  vector[N_tr] y_mat_tr[T]; // response, cover ts
  vector[N_tr] y_preb; // pre-burn cover
  int NN_tr[N_tr,4];
  int NNind_tr[N_tr,4];
  real NN01_tr[N_tr,4];
  /*
  int N_nb;   // total number of neighbors (for indexing)
  int nbInd[N_nb]; // indices for <vector> u
  int fInd[N_nb]; // index for spatial functions [1:5], indexing for neighbor positions
  int numNb[N];  // number of neighbors for each pixel (for segment())
  int pos[N]; // position vector (for segement())
  */
} transformed data {
  real dx = 5; 
  real dy = 5;
  //matrix[T,N_tr] y = y_mat_tr';
} parameters {
  real alpha;
  real<lower=0> r; // pixel-level growth/effect of pixel on itself
  // real<lower=0> r_mu;
  // real<lower=0> r_sigma;
  real<lower=0> k;
  // real<lower=0> k_mu;
  // real k_phi;
  // real<lower=0> k_sigma;
  vector<lower=0>[N_tr] D; // dispersion parameter
  vector<lower=0>[T] eta;  // observation error
  vector<lower=0>[T-1] gamma;  // observation error
  vector<lower=0>[N_tr] z_init;    // initial state vector
  vector<lower=0>[N] z_pred;    // initial state vector
} transformed parameters {
  vector[N_tr] Z[T]; // latent states
  matrix[N_tr,N_tr] H_tr;
  H_tr = H_mat(D, r, dx, dy, NN_tr, NN01_tr, N_tr);

  Z[1] =z_init;
  
    /*
  vector[rows(csr_extract_w(H'))] w; 
  int v[4060]; 
  int u[N+1]; 
  w = csr_extract_w(H');
  v = csr_extract_v(H');
  u = csr_extract_u(H');
  */
  for(t in 2:T)
      Z[t] = Z[t-1]+(H_tr*Z[t-1])+r[1]*Z[t-1] - r*(Z[t-1] .* Z[t-1])/k;// + gamma[t-1];
      //(csr_matrix_times_vector(N, N, w, v, u, Z[t-1]'))' 
  
} model {
  //priors
  alpha ~ normal(0,1);
  
  // r_mu ~ normal(0,.1);
  // r_sigma ~ exponential(5);
  r ~ normal(0,5);
  // 
  // k_mu ~ normal(0,5);
  // k_phi ~ normal(0,1);
  // k_sigma ~ exponential(10);

  //D ~ inv_gamma(10,20);
  D ~ normal(0,5);
  
  eta ~ exponential(1); 
  gamma ~ normal(0,.1);
  z_init ~ normal(0,5);
  //likelihoods
  //k ~ normal(k_mu*y_preb, k_sigma);
  y_mat_tr[1] ~ normal(z_init,eta[1]);//normal(z_init,1/(z_init*eta[1]));
  y_mat[1] ~ normal(z_pred,eta[1]); 
   for(t in 2:T){
      y_mat_tr[t] ~ normal(alpha+Z[t],eta[t]);//normal(Z[t],1/(Z[t]*eta[t]));
   }
} generated quantities {
  vector<lower=0>[N] y_pred[T];
  matrix[N,N] H;
  vector[N] D_mean = rep_vector(mean(D),N);
  H = H_mat(D_mean, r[1], dx, dy, NN, NN01, N);
  for(i in 1:N)
    y_pred[1,i] = z_pred[i];//fabs(normal_rng(z_init[i],eta[1]));
  for(t in 2:T){
    //for(i in 1:N)
      y_pred[t] = y_pred[t-1]+(H*y_pred[t-1])+r[1]*y_pred[t-1] - r[1]/r[2]*(y_pred[t-1] .* y_pred[t-1]);//+ gamma[t-1];
    }
}
