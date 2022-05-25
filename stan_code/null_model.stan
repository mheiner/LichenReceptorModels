data {
  int<lower=0> L; // number of elements
  int<lower=0> n; // number of samples
  int<lower=0> K; //number of profiles
  matrix[L,K] alpha; // prior alphas for baseline profile means/sds
  matrix[L,K] beta; // prior betas for baseline profile means/sds
  vector[L] a;         // prior alpha for cv
  vector[L] b;         // prior beta for cv
  real flim;     // upper bound for uniform draws for f matrix
  
  // Inserting bdl/missing value pieces
  int n_obs;
  int n_mis;
  int n_bdl;
  int indx_obs[n_obs]; // indexes a flattened (column major order) L by n matrix
  int indx_mis[n_mis]; // indexes a flattened (column major order) L by n matrix
  int indx_bdl[n_bdl]; // indexes a flattened (column major order) L by n matrix
  vector<lower=0>[n_obs] y_obs;
  vector<lower=0>[n_bdl] dl_censor_vals;
}

transformed data {
  vector[n_obs] ly_obs = log(y_obs);
  vector[n_bdl] ldl_censor_vals = log(dl_censor_vals);
  matrix[L,K] logbeta = log(beta);
}

parameters {
  matrix[L-1, K] x; // logit of lam
  row_vector<upper=log(flim)>[n] lF0;
  vector<lower=0>[L] cv;
  
  vector<lower=0>[n_bdl] logamount_bdl;
  vector[n_mis] ly_mis;
}

transformed parameters {
  vector[n*L] ly; // complete data vector

  ly[indx_obs] = ly_obs;
  ly[indx_mis] = ly_mis;
  ly[indx_bdl] = (ldl_censor_vals - logamount_bdl);

  // llam_vec[indx_non0Lam] = lx;
  // llam_vec[indx_0Lam] = rep_vector(negative_infinity(), n0Lam);
}

model {
  matrix[L,n] lLF;
  vector[L*n] sig;


  vector[K] lxdenom;
  matrix[L,K] llam = rep_matrix(negative_infinity(), L, K);
  for(k in 1:K) {
    lxdenom[k] = log1p_exp(log_sum_exp(x[,k]));
    for(ell in 1:(L-1)) {
      llam[ell, k] = x[ell,k] - lxdenom[k];
    }
    llam[L, k] = -lxdenom[k];
  }

  for(k in 1:K) { // Generalized Dirichlet Prior of Lingwall et. al. (2008)
    target += -sum(alpha[, k])*log_sum_exp(llam[, k] + logbeta[, k]); // beta is rate parameter
    target += sum(alpha[, k] .* llam[, k]); // Jacobian for log(lam) built in here
  }

  exp(lF0) ~ uniform(0, flim);  // baseline contributions, row vector
  target += sum(lF0);
  
  cv ~ gamma(a, b);

  for (j in 1:L){
    for(i in 1:n){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF0[i]);
    }
  }
  
  sig = sqrt(log(square(to_vector(rep_matrix(cv, n))) + 1.0));
  
  ly ~ normal( to_vector(lLF), sig );
}

generated quantities { // for monitoring
  // matrix[L,K] llam = to_matrix(llam_vec, L, K);
  matrix[L,K] lam;
  matrix[L,n] lLF;
  vector[L*n] sig;
  vector[n_obs] llik;
  
  vector[K] lxdenom;
  matrix[L,K] llam = rep_matrix(negative_infinity(), L, K);
  for(k in 1:K) {
    lxdenom[k] = log1p_exp(log_sum_exp(x[,k]));
    for(ell in 1:(L-1)) {
      llam[ell, k] = x[ell,k] - lxdenom[k];
    }
    llam[L, k] = -lxdenom[k];
  }

  lam = exp(llam);

  for (i in 1:n){
    for(j in 1:L){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF0[i]);
    }
  }
  
  sig = sqrt(log(square(to_vector(rep_matrix(cv, n))) + 1.0));
  for (i in 1:n_obs) {
    llik[i] = normal_lpdf(ly_obs[i] | to_vector(lLF)[indx_obs[i]], sig[indx_obs[i]]);
  }
}