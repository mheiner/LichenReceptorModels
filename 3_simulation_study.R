rm(list=ls())

ARGS = commandArgs(trailingOnly = TRUE) # command-line arguments
# ARGS = c(3, 0, 1, 1.0, 50, 1234, "sparseT", "Id2")

K0 = as.numeric(ARGS[1]) # Number of original sources for fitting the model (number of source priors -- 3 or 4 including baseline)
Ke = as.numeric(ARGS[2]) # number of "extra" sources 0, 1 (natural), 2 (natural, anthro)
src_present = as.numeric(ARGS[3]) # Number of actual sources (beyond baseline) actually present in the data -- 1 or 3
inflation_factor = as.numeric(ARGS[4]) # 0.2, 1, 3
sparsity = as.numeric(ARGS[5]) # Level of sparsity in the source contributions -- 0% or 50%
seed = as.numeric(ARGS[6])
model = as.character(ARGS[7]) # one of "base" "sparse" "baseT" "sparseT" "RCsparseT" or "sparseDL"
idf_rule = as.character(ARGS[8]) # one of "original" and "Id2"

library("parallel")
library("rstan")

N = 100
effective_zero = 1e-12

set.seed(seed)

min_alpha = 2.0
source('helperScripts/all_profiles_nitrogen_modified_specificInflation_unimodal.R')

el_use = 1:ncol(alpha_matrix) # Use all 25 elements
(L <- length(el_use))

# Create the Lambda matrix that will be used to generate the data (make it match the full fitting Lambda)
if (src_present == 1) {
  
  src_use <- c(1, 6) # Use baseline and 1 other source (in this case unpaved road dust)
  x <- t(alpha_matrix[src_use, el_use] / beta_matrix[src_use, el_use])
  Lam_truth = x / rep(colSums(x), each = nrow(x))
  colSums(Lam_truth)
  
  # Fix zeroes in Lam_truth to satisfy identifiability conditions
  for (k in 1:ncol(Lam_truth) ) {
    ord = order(Lam_truth[,k])
    Lam_truth[ord[1:(ncol(Lam_truth)-1)],k] = effective_zero
  }
  
} else if (src_present == 3) {
  
  src_use <- c(1, 2, 4, 6) # Use baseline and 4 other sources (in this case playa, brake, unpaved dust)
  x <- t(alpha_matrix[src_use, el_use] / beta_matrix[src_use, el_use])
  Lam_truth = x / rep(colSums(x), each = nrow(x))
  colSums(Lam_truth)
  
  # Fix zeroes in Lam_truth to satisfy identifiability conditions
  for (k in 3:4 ) {
    ord = order(Lam_truth[,k])
    Lam_truth[ord[1:(ncol(Lam_truth)-1)],k] = effective_zero
  }
  Lam_truth[c(8,11,19),1] = effective_zero
  Lam_truth[c(12,16,17),2] = effective_zero
}

# Re-normalize
(Lam_truth = Lam_truth / rep(colSums(Lam_truth), each = nrow(Lam_truth)))
colSums(Lam_truth)

source("helperScripts/test_profiles.R")
rankcheck(Lam_truth, tol=100*effective_zero)

# Separate out the baseline profile and the non-baseline profile(s)
Lam_0 = Lam_truth[, 1]
Lam = Lam_truth[, -1]

# Generate the source contribution matrices f_0 and F_mat
flim <- 120000 # PPM
f_0 <- runif(N, 50000, flim) # PPM

F_mat <- matrix(rgamma(src_present*N, 2, scale = 10000), nrow=src_present, ncol=N) # all in PPM


# Replace 0 or 50% with zeroes at random
if(sparsity > 0) {
  forced_zeroes <- sample(1:length(F_mat), floor(0.01*sparsity*length(F_mat)), replace = F)
  F_mat[forced_zeroes] <- effective_zero
} else {
  F_mat <- F_mat
}

if (grepl("T", model)) {
  cv = runif(L, 0.1, 0.8)
  (sig = sqrt(log(cv^2 + 1.0)))
  Y <- matrix( exp( rep(sig, times=N) * rt(L*N, df=4.5) + log((Lam_0 %*% t(f_0)) + (Lam %*% F_mat)) ), nrow=L)
} else {
  (cv = runif(L, 0.2, 1.2))
  (sig = sqrt(log(cv^2 + 1.0)))
  Y <- matrix( exp(rnorm( L*N, log((Lam_0 %*% t(f_0)) + (Lam %*% F_mat)), rep(sig, times=N))), nrow=L)
}

# max(abs(Lam_truth %*% rbind(f_0, F_mat) - ((Lam_0 %*% t(f_0)) + (Lam %*% F_mat))))

dim(Y)

### Randomly assign missing values and determine BDL values (we will consider values <= 5th percentile for that element to be BDL)
n_mis = 150
Y_new = Y
indx_mis = sort(sample(1:(dim(Y)[1]*dim(Y)[2]), n_mis, replace=FALSE))
Y_new[indx_mis] = NA

detec_limits = apply(Y_new, 1, FUN = function(x) quantile(probs=0.05, x, na.rm=TRUE))
detec_limits
bdl = (Y_new <= matrix(rep(detec_limits, times=ncol(Y_new)), ncol=ncol(Y_new)))
Y_new[bdl] = NA
indx_bdl = which(bdl)
indx_obs = which(!is.na(Y_new))

y_new = as.vector(Y_new)
y_obs = Y_new[!is.na(y_new)]
(n_obs = length(y_obs))
(n_bdl = sum(bdl, na.rm=TRUE))

all_detec_limits = rep(detec_limits, times=N)
dl_censor_vals = all_detec_limits[indx_bdl]




# Use either 6 or 7 profiles as priors for the model
if (K0 == 4) {
  src_use <- c(1, 2, 4, 6) # Baseline, playa, brake, unpaved.
} else if (K0 == 3) {
  # Here, we exclude "unpaved" from the priors, even though we know it is present in generating the data
  src_use <- c(1, 2, 4) # Baseline, playa, brake, (unpaved excluded).
}

### Prior
if (Ke == 0) {
  Alpha = t(alpha_matrix[src_use, el_use])
  Beta = t(beta_matrix[src_use, el_use])
} else if (Ke == 1) {
  Alpha = t(rbind(alpha_matrix[src_use, el_use], empty_alpha))
  Beta = t(rbind(beta_matrix[src_use, el_use], natural_beta))
} else if (Ke == 2) {
  Alpha = t(rbind(alpha_matrix[src_use, el_use], empty_alpha, empty_alpha))
  Beta = t(rbind(beta_matrix[src_use, el_use], natural_beta, anthropogenic_beta))
}

(K = K0 + Ke)

dim(Alpha); dim(Beta)
(src_names = gsub("_beta", "", colnames(Beta)))

## model identifiability zeros
x = Alpha / Beta
(Lam_prior = x / rep(colSums(x), each=nrow(x)))
colnames(Lam_prior) = src_names
colSums(Lam_prior)

if (idf_rule == "original") {
  
  ## original identifiability conditions
  for (k in 3:K0) {
    ord = order(Lam_prior[,k])
    Lam_prior[ord[1:(K-1)],k] = 0.0
  }
  if (K0 == 3) {
    Lam_prior[c(8,11),1] = 0.0
    Lam_prior[c(12,16),2] = 0.0
    # Natural is profile K-1, Anthro is K:
    if (Ke == 1) {
      Lam_prior[c(8,11,19),1] = 0.0
      Lam_prior[c(12,16,17),2] = 0.0
      Lam_prior[c(13,14,19), K0+1] = 0.0 # Natural
    }
    if (Ke == 2) {
      Lam_prior[c(8,11,19,17),1] = 0.0
      Lam_prior[c(12,16,17,25),2] = 0.0
      Lam_prior[c(13,14,19,20), K0+1] = 0.0 # Natural
      Lam_prior[order(Lam_prior[,K])[1:(K-1)], K0+2] = 0.0 # Anthropogenic
    }
  }
  if (K0 == 4) {
    Lam_prior[c(8,11,19),1] = 0.0
    Lam_prior[c(12,16,17),2] = 0.0
    # Natural is profile K-1, Anthro is K:
    if (Ke == 1) {
      Lam_prior[c(8,11,19,17),1] = 0.0
      Lam_prior[c(12,16,17,25),2] = 0.0
      Lam_prior[c(13,14,19,20), K0+1] = 0.0 # Natural
    }
    if (Ke == 2) {
      Lam_prior[c(8,11,19,17,21),1] = 0.0
      Lam_prior[c(12,16,17,25,21),2] = 0.0
      Lam_prior[c(13,14,19,20,21), K0+1] = 0.0 # Natural
      Lam_prior[order(Lam_prior[,K])[1:(K-1)], K0+2] = 0.0 # Anthropogenic
    }
  }
  
} else if (idf_rule == "Id2") {
  indx_zero = read.csv(file="data/hardZeroId2.csv", header=TRUE)
  nz = K - 1
  for (k in 1:ncol(Lam_prior)) {
    Lam_prior[indx_zero[1:nz, grep( src_names[k], colnames(indx_zero))], k] = 0
  }
}

rankcheck(Lam_prior, 0)

dim(Lam_prior)
indx_non0Lam = which(Lam_prior > 100*effective_zero) # vector indexing of the L by K Lambda matrix
indx_0Lam = which(Lam_prior <= 100*effective_zero)
(n_non0Lam = length(indx_non0Lam))
(n0Lam = length(indx_0Lam))
all.equal(sort(c(indx_non0Lam, indx_0Lam)), 1:(L*K))

(n0Lam_each = sum(Lam_prior[,1] <= 100*effective_zero))
stopifnot( all(apply(Lam_prior, 2, function(x) sum(x <= 100*effective_zero)) == n0Lam_each) )
non0Lam_indx = apply(Lam_prior, 2, function(x) which(x > 100*effective_zero))


## Other prior settings
if (model %in% c("base", "sparse")) {
  mean_cv = rep(0.4, L)
  sd_cv = rep(0.15, L)
  a = mean_cv^2 / sd_cv^2
  b = mean_cv / sd_cv^2
} else if (model %in% c("baseT", "sparseT")) {
  a = rep(1.2, L)
  b = rep(4.8, L) # "cv" on log-t with low df is very high, pretty close to 1-1 with scale of t-likelihood
}

a / b
sqrt(a / b^2)


gamma_0 <- 5 * 10e3 # PPM


## data
dat_stan = list(y_obs = y_obs, L = L, n = N, K = K, Ke = Ke, alpha = Alpha, beta = Beta, 
				a = a, b = b, flim = flim,
                gamma_0 = gamma_0, n_obs = n_obs, n_mis = n_mis, n_bdl = n_bdl,
                indx_obs = indx_obs, indx_mis = indx_mis, indx_bdl = indx_bdl, 
                dl_censor_vals = dl_censor_vals,
                n0Lam_each = n0Lam_each, non0Lam_indx = non0Lam_indx)

if (model == "base") {
	stanfile = 'stan_code/base_model.stan'
} else if (model == "baseT") {
    stanfile = 'stan_code/baseT_model.stan'
} else if (model == "sparse") {
    stanfile = 'stan_code/sparse_model_margF.stan'
} else if (model == "sparseT") {
    stanfile = 'stan_code/sparseT_model_margF.stan'
}


n_chains =  3
(n_cores = min( parallel::detectCores(logical=FALSE), n_chains ))

# Sample from the Posterior Distribution

fit = rstan::stan(
  file = stanfile,
  data = dat_stan,
  chains = n_chains,
  cores = n_cores,
  # warmup = 2000, iter = 3000,
  warmup = 100, iter = 200, # for testing
  control = list(max_treedepth = 15,
                 adapt_delta = 0.95),
  pars = c("x", "lx", "ly", "ly_mis", "logamount_bdl", "llam", "llam_vec", 
           "llam0", "phi", "lF", "lF0", "lF1", "lFe", "lxdenom", "sig"), # "lLF", "llik"
  seed = seed,
  include = FALSE
)

summ = summary(fit)

# Export summaries
median_contributions <- summary(fit, pars='F')$summary[ , "50%"]
mean_profiles <- summary(fit, pars='lam')$summary[,'mean']

dim(mean_profiles) <- c(K, L)
dim(median_contributions) <- c(N, K)

median_contributions <- as.data.frame(median_contributions)

median_lFL = summary(fit, pars="lLF")$summary[,"50%"] # arranged with sample index changing first
dim(median_lFL) = c(N, L)

rhats <- summary(fit)$summary[,"Rhat"]
head(rhats[order(rhats, decreasing = TRUE)], n=30)

cvs <- summary(fit, pars='cv')$summary[,'mean']
cv_true = cv

chains = extract(fit, pars=c("lLF"), include=FALSE)


# Save model fits and results
save(chains, 
     file = paste0('postsim/', model, '_sim_chains', '_idf', idf_rule, '_K0', K0, '_Ke', Ke, '_', src_present, 'present_', inflation_factor, 'infl_', sparsity, 'sparse_seed', seed, '.rda'))

# save(fit, 
#      file = paste0('postsim/', model, '_sim_fit', '_idf', idf_rule, '_K0', K0, '_Ke', Ke, '_', src_present, 'present_', inflation_factor, 'infl_', sparsity, 'sparse_seed', seed, '.rda'))

save(rhats, summ, median_lFL, median_contributions, mean_profiles, 
     cvs, cv_true,
     f_0, F_mat, Lam_truth, dat_stan, Y, Y_new, Lam_prior,
     file = paste0('postsim/', model, '_sim_results', '_idf', idf_rule, '_K0', K0, '_Ke', Ke, '_', src_present, 'present_', inflation_factor, 'infl_', sparsity, 'sparse_seed', seed, '.rda'))

quit(save="no")