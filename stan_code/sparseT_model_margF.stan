data {
  int<lower=0> L; // number of elements
  int<lower=0> n; // number of samples ( n*L = n_obs + n_mis + n_bdl )
  int<lower=0> K; //number of profiles
  matrix[L,K] alpha; // prior alphas for baseline profile means/sds
  matrix[L,K] beta; // prior betas for baseline profile means/sds
  vector[L] a;         // prior alpha for cv (beta prior)
  vector[L] b;         // prior beta for cv (beta prior)
  real flim;     // upper bound for uniform draws for f matrix
  real gamma_0;
  // Inserting bdl/missing value pieces
  int n_obs;
  int n_mis;
  int n_bdl;
  int indx_obs[n_obs]; // indexes a flattened (column major order) L by n matrix
  int indx_mis[n_mis]; // indexes a flattened (column major order) L by n matrix
  int indx_bdl[n_bdl]; // indexes a flattened (column major order) L by n matrix
  vector<lower=0>[n_obs] y_obs;
  vector<lower=0>[n_bdl] dl_censor_vals;
  
  int<lower=0> n0Lam_each;
  int<lower=0> non0Lam_indx[L-n0Lam_each, K];
}

transformed data {
  vector[n_obs] ly_obs = log(y_obs);
  vector[n_bdl] ldl_censor_vals = log(dl_censor_vals);
  matrix[L,K] logbeta = log(beta);
}

parameters {
  matrix[L-n0Lam_each-1, K] x; // logit of lam

  row_vector<upper=log(flim)>[n] lF0;
  matrix[K-1,n] lF1;
  vector<lower=0, upper=1>[L] cv;
  real<lower=4> nu;

  row_vector<lower=0>[n] gamma;
  vector<lower=0>[n_bdl] logamount_bdl;
  vector[n_mis] ly_mis;
}

transformed parameters {
  vector[n*L] ly; // complete data vector

  ly[indx_obs] = ly_obs;
  ly[indx_mis] = ly_mis;
  ly[indx_bdl] = (ldl_censor_vals - logamount_bdl);
}

model {
  matrix[K,n] lF;
  matrix[L,n] lLF;
  vector[L*n] sig;
  
  vector[K] lxdenom;
  matrix[L,K] llam = rep_matrix(negative_infinity(), L, K);
  for(k in 1:K) {
    lxdenom[k] = log1p_exp(log_sum_exp(x[,k]));
    for(ell in 1:(L-n0Lam_each-1)) {
      llam[non0Lam_indx[ell, k], k] = x[ell,k] - lxdenom[k];
    }
    llam[non0Lam_indx[L-n0Lam_each, k], k] = -lxdenom[k];
  }

  for(k in 1:K) { // Generalized Dirichlet Prior of Lingwall et. al. (2008)
    target += -sum(alpha[non0Lam_indx[,k], k])*log_sum_exp(llam[non0Lam_indx[,k], k] + logbeta[non0Lam_indx[,k], k]); // beta is rate param
    target += sum(alpha[non0Lam_indx[,k], k] .* llam[non0Lam_indx[,k], k]); // Jacobian for x built in here
  }

  gamma ~ exponential(1 / (2 * gamma_0));

  exp(lF0) ~ uniform(0, flim);  // baseline contributions, row vector
  target += sum(lF0);
  
  // to_vector(exp(lF1)) ~ gamma(1.0/(K-1), to_vector(rep_matrix(inv(gamma), K-1))); // everything excluding baseline
  to_vector(exp(lF1)) ~ gamma(0.5, to_vector(rep_matrix(inv(gamma), K-1))); // everything excluding baseline
  target += sum(to_vector(lF1));

  cv ~ beta(a, b);
  nu ~ inv_gamma(19, 200);

  // transformations
  lF = append_row(lF0, lF1);
  for (i in 1:n){
    for(j in 1:L){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF[,i]);
    }
  }

  sig = sqrt(log(square(to_vector(rep_matrix(cv, n))) + 1.0)); 
  
  // data
  ly ~ student_t(nu, to_vector(lLF), sig);
}

generated quantities { // for monitoring
  matrix[L,K] lam;
  matrix[K,n] F;
  matrix[K,n] lF;
  matrix[L,n] lLF;
  vector[L*n] sig;
  vector[n_obs] llik;

  vector[K] lxdenom;
  matrix[L,K] llam = rep_matrix(negative_infinity(), L, K);
  for(k in 1:K) {
    lxdenom[k] = log1p_exp(log_sum_exp(x[,k]));
    for(ell in 1:(L-n0Lam_each-1)) {
      llam[non0Lam_indx[ell, k], k] = x[ell,k] - lxdenom[k];
    }
    llam[non0Lam_indx[L-n0Lam_each, k], k] = -lxdenom[k];
  }

  lam = exp(llam);
  lF = append_row(lF0, lF1);
  F = exp(lF);
  
  for (i in 1:n){
    for(j in 1:L){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF[,i]);
    }
  }
  
  sig = sqrt(log(square(to_vector(rep_matrix(cv, n))) + 1.0));
  for (i in 1:n_obs) {
    llik[i] = student_t_lpdf(ly_obs[i] | nu, to_vector(lLF)[indx_obs[i]], sig[indx_obs[i]]);
  }
}