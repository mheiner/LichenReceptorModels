
data {
  int<lower=0> L; // number of elements
  int<lower=0> n; // number of samples
  int<lower=0> K; //number of profiles
//  matrix[L,n] Y; //observed data
  matrix[L,K] alpha; // prior alphas for baseline profile means/sds
  matrix[L,K] beta; // prior betas for baseline profile means/sds
  vector[L] a;         // prior alpha for cv
  vector[L] b;         // prior beta for cv
  real flim;     // upper bound for uniform draws for f matrix
  
  // Inserting bdl/missing value pieces
  int n_obs;
  int n_mis;
  int n_bdl;
  int indx_obs[n_obs]; // indexes a flattened (column major order) n by L matrix
  int indx_mis[n_mis]; // indexes a flattened (column major order) n by L matrix
  int indx_bdl[n_bdl]; // indexes a flattened (column major order) n by L matrix
  vector<lower=0>[n_obs] y_obs;
  vector<lower=0>[n_bdl] dl_censor_vals;
  
  int<lower=0> n0Lam; // number of hard zeros in profiles
  int<lower=0> n_non0Lam; // L*K - n0Lam
  int indx_0Lam[n0Lam];
  int indx_non0Lam[n_non0Lam];
}

transformed data {
  vector[n_obs] ly_obs = log(y_obs);
  vector[n_bdl] ldl_censor_vals = log(dl_censor_vals);
  vector[n_non0Lam] alpha_vec = to_vector(alpha)[indx_non0Lam];
  vector[n_non0Lam] beta_vec = to_vector(beta)[indx_non0Lam];
}

parameters {
  vector[n_non0Lam] lx;
  row_vector<upper=log(flim)>[n] lF0;
  matrix[K-1,n] lF1;
  vector<lower=0>[L] cv;
  
  vector<lower=0>[n_bdl] logamount_bdl;
  vector[n_mis] ly_mis;
}

transformed parameters {
  vector[n*L] ly; // complete data vector
  vector[L*K] llam_vec;

  ly[indx_obs] = ly_obs;
  ly[indx_mis] = ly_mis;
  ly[indx_bdl] = (ldl_censor_vals - logamount_bdl);

  llam_vec[indx_non0Lam] = lx;
  llam_vec[indx_0Lam] = rep_vector(negative_infinity(), n0Lam);
}

model {
  matrix[K,n] lF;
  matrix[L,n] lLF;
  vector[L*n] sig;

  matrix[L,K] llam = to_matrix(llam_vec, L, K);
  for(k in 1:K) {
    llam[,k] = llam[,k] - log_sum_exp(llam[,k]);
  }

  exp(lx) ~ gamma(alpha_vec, beta_vec);
  target+= sum(lx);

  exp(lF0) ~ uniform(0, flim);  // baseline contributions, row vector
  target += sum(lF0);
  
  to_vector(lF1) ~ normal(log(1), sqrt(log(square(3) + 1.0))); // non-baseline

  cv ~ gamma(a, b);

  lF = append_row(lF0, lF1);
  for (j in 1:L){
    for(i in 1:n){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF[,i]);
    }
  }
  
  sig = sqrt(log(square(to_vector(rep_matrix(cv, n))) + 1.0));
  
  ly ~ normal( to_vector(lLF), sig );
}

generated quantities { // for monitoring
  matrix[L,K] llam = to_matrix(llam_vec, L, K);
  for(k in 1:K) {
    llam[,k] = llam[,k] - log_sum_exp(llam[,k]);
  }
  matrix[L,K] lam = exp(llam);
  matrix[K,n] F = exp(append_row(lF0, lF1));
  matrix[K,n] lF = log(F);

  matrix[L,n] lLF;
  for (j in 1:L){
    for(i in 1:n){
      lLF[j,i] = log_sum_exp(to_vector(llam[j,]) + lF[,i]);
    }
  }
}