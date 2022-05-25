ARGS = commandArgs(trailingOnly = TRUE)
mod = as.character(ARGS[1])
idf_rule = as.character(ARGS[2])
K0 = as.numeric(ARGS[3])
Ke = as.numeric(ARGS[4])
inflation_factor = as.numeric(ARGS[5]) # replaces the assignment at beginning of "helperScripts/all_prof...R" (be sure that is commented out)
seed = as.numeric(ARGS[6])


## model settings used in paper
# mod = "sparseT"
# idf_rule = "Id2"
# K0 = 5
# Ke = 2
# inflation_factor = 1.0
# seed = 42702


library("parallel")
library("rstan")


# Load in the raw rhizoplaca lichen data
load(file = 'data/dat_rhizICP.rda')

min_alpha = 5.0 # unimodal gammas in profile priors
source('helperScripts/all_profiles_nitrogen_modified_specificInflation_unimodal.R')

L = dat_stan$L

## prior on error coefficient of variation
if (mod %in% c("base", "sparse", "null")) {
  mean_cv = rep(0.4, L)
  sd_cv = rep(0.15, L)
  a = mean_cv^2 / sd_cv^2
  b = mean_cv / sd_cv^2
} else if (mod %in% c("baseT", "sparseT", "nullT")) {
  a = rep(1.2, L)
  b = rep(4.8, L) # "cv" on log-t with low df is very high, pretty close to 1-1 with scale of t-likelihood
}

## use baseline, playa, brake, exhaust, and unpaved profiles.
if (K0 == 5) {
  src_use = c(1,2,4,5,6)
} else if (K0 == 1) {
  src_use = 1 # baseline only for null model
}

if (Ke == 0) {
  Alpha = rbind(alpha_matrix[src_use,, drop=FALSE])
  Beta = rbind(beta_matrix[src_use,, drop=FALSE])
} else if (Ke == 1) {
  Alpha = rbind(alpha_matrix[src_use,, drop=FALSE], empty_alpha)
  Beta = rbind(beta_matrix[src_use,, drop=FALSE], natural_beta)
} else if (Ke == 2) {
  Alpha = rbind(alpha_matrix[src_use,, drop=FALSE], empty_alpha, empty_alpha)
  Beta = rbind(beta_matrix[src_use,, drop=FALSE], natural_beta, anthropogenic_beta)
}

(K = K0 + Ke) # total number of sources, including baseline and Ke extras

el_use = 1:ncol(Alpha)
(L = length(el_use))
x = t(Alpha[, el_use, drop=FALSE] / Beta[, el_use, drop=FALSE])
(Lam_prior = x / rep(colSums(x), each=nrow(x)))

(src_names = gsub("_beta", "", rownames(Beta)))
colnames(Lam_prior) = src_names

## Identifiability rule (structured zeros in profile matrix)
if (idf_rule == "original" & (K > 1)) {
  
  ## original identifiability conditions
  for (k in 3:K0) {
    ord = order(Lam_prior[,k])
    Lam_prior[ord[1:(K-1)],k] = 0.0
  }
  Lam_prior[c(8,11,12,17),1] = 0.0 # add in element 21 if more than 5 profiles
  Lam_prior[c(12,16,17,25),2] = 0.0 # add in element 21 if more than 5 profiles
  if (K0 == 5 & Ke == 1) {
    Lam_prior[21,1] = 0.0
    Lam_prior[21,2] = 0.0
    Lam_prior[c(13,14,19,20,21), 6] = 0.0
  }
  if (K0 == 5 & Ke == 2) {
    Lam_prior[c(8,11,12,17,19,21),1] = 0.0 # Baseline
    Lam_prior[c(12,16,17,19,21,25),2] = 0.0 # Playa
    Lam_prior[c(13,14,19,20,21,25),6] = 0.0 # Natural
    Lam_prior[order(Lam_prior[,7])[1:(K-1)],7] = 0.0 # Anthropogenic
  }
  
} else if (idf_rule == "Id2" & (K > 1)) {
  
  indx_zero = read.csv(file="data/hardZeroId2.csv", header=TRUE)
  nz = K - 1
  for (k in 1:ncol(Lam_prior)) {
    Lam_prior[indx_zero[1:nz, grep( src_names[k], colnames(indx_zero))], k] = 0
  }

}

Lam_prior
(Lam_prior = Lam_prior / rep(colSums(Lam_prior), each=nrow(Lam_prior)))
colSums(Lam_prior)

source("helperScripts/test_profiles.R")
rankcheck(Lam_prior, 0)


#############################################################################################################

# Make sure prior values are all assigned properly for the model:
dat_stan$alpha = t(Alpha)
dat_stan$beta = t(Beta)

dat_stan$a = a
dat_stan$b = b
dat_stan$K = K

dat_stan$Lam_prior = Lam_prior

if (grepl("null", mod)) {
  dat_stan$flim = 500e3
}

dat_stan$Ke = Ke

dat_stan$gamma_0 = 3.0 * 10e3 # in PPM

dat_stan$indx_non0Lam = which(Lam_prior > 0.0) # vector indexing of the L by K Lambda matrix
dat_stan$indx_0Lam = which(Lam_prior == 0.0)
dat_stan$n_non0Lam = length(dat_stan$indx_non0Lam)
dat_stan$n0Lam = length(dat_stan$indx_0Lam)
all.equal(sort(c(dat_stan$indx_non0Lam, dat_stan$indx_0Lam)), 1:(dat_stan$L*dat_stan$K))

(dat_stan$n0Lam_each = sum(Lam_prior[,1] == 0))
stopifnot( all(apply(Lam_prior, 2, function(x) sum(x == 0)) == dat_stan$n0Lam_each) )
dat_stan$non0Lam_indx = apply(dat_stan$Lam_prior, 2, function(x) which(x > 0))




n_chains = 4
# n_chains = 3 # for sensitivity analysis
(n_cores = min( parallel::detectCores(logical=FALSE), n_chains ))


if (mod == "base") {
  file1 = "stan_code/base_model.stan"
} else if (mod == "baseT") {
  file1 = "stan_code/baseT_model.stan"  
} else if (mod == "sparse") {
  file1 = "stan_code/sparse_model_margF.stan"
} else if (mod == "sparseT") {
  file1 = "stan_code/sparseT_model_margF.stan"
} else if (mod == "null") {
  file1 = "stan_code/null_model.stan"
} else if (mod == "nullT") {
  file1 = "stan_code/nullT_model.stan"
}


fit = rstan::stan(
  file = file1, data = dat_stan,
  chains = n_chains, cores = n_cores,
  # iter = 15000, warmup = 10000,
  iter = 200, warmup = 100, # for testing
  control = list(max_treedepth = 15,
                 adapt_delta = 0.95),
  pars = c("x", "lx", "ly", "ly_mis", "logamount_bdl", 
           "llam", "llam_vec", "llam0", "phi", 
           "lF", "lF0", "lF1", "lFe", "lxdenom", 
           "sig"), #, "llik", "lLF"),
  include=FALSE,
  seed = seed
)

summ = summary(fit)
rhats = summ$summary[,"Rhat"]

chains = extract(fit, pars=c("lLF"), include=FALSE)


# Save model fits and results
save(chains, rhats,
     file = paste0('postsim/chains_', mod, 'Model_idf_', idf_rule, '_K0', K0, '_Ke', Ke, '_inf', inflation_factor, "_seed", seed, '.rda'))

save(fit, dat_stan, 
     file=paste0('postsim/fit_', mod, 'Model_idf_', idf_rule, '_K0', K0, '_Ke', Ke, '_inf', inflation_factor, "_seed", seed, '.rda'))

save(summ, rhats, dat_stan, 
     file=paste0('postsim/results_', mod, 'Model_idf_', idf_rule, '_K0', K0, '_Ke', Ke, '_inf', inflation_factor, "_seed", seed, '.rda'))

quit(save="no")