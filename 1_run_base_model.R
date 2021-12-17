
ARGS = commandArgs(trailingOnly = TRUE)
idf_rule = ARGS[1]
K = as.numeric(ARGS[2])
inflation_factor = as.numeric(ARGS[3]) # replaces the assignment at beginning of "all_prof...R" (be sure to comment that out)
seed = as.numeric(ARGS[4])

## sensitivity settings
# idf_rule = "original"
# K = 5
# inflation_factor = 1.0 # replaces the assignment at beginning of "all_prof...R"
# seed = 1

# Load in the raw rhizoplaca lichen data
load(file = 'data/dat_rhizICP.rda')

# Run the file below to access prior information for each profile
source('helperScripts/all_profiles_nitrogen_modified_specificInflation.R')

L = dat_stan$L

## prior on error coefficient of variation
mean_cv = rep(0.5, L)
sd_cv = rep(0.2, L)

a = mean_cv^2 / sd_cv^2
b = mean_cv / sd_cv^2


# Use baseline, playa, brake, exhaust, and unpaved profiles.
# Impose zeroes in specific places in order to ensure identifiability conditions are satisfied.
src_use = c(1,2,4,5,6)

if (K == 5) {
  alpha_matrix = rbind(alpha_matrix[src_use,])
  beta_matrix = rbind(beta_matrix[src_use,])
} else if (K == 6) {
  alpha_matrix = rbind(alpha_matrix[src_use,], natural)
  beta_matrix = rbind(beta_matrix[src_use,], empty_beta)
} else if (K == 7) {
  alpha_matrix = rbind(alpha_matrix[src_use,], natural, anthropogenic)
  beta_matrix = rbind(beta_matrix[src_use,], empty_beta, empty_beta)
}


el_use = 1:ncol(alpha_matrix)
# (K = nrow(alpha_matrix))
(L = length(el_use))
x = t(alpha_matrix[, el_use] / beta_matrix[, el_use])
(Lam_all = x / rep(colSums(x), each=nrow(x)))


if (idf_rule == "original") {
  
  ## original identifiability conditions
  for (k in 3:5) {
    ord = order(Lam_all[,k])
    Lam_all[ord[1:(K-1)],k] = 0.0
  }
  Lam_all[c(8,11,12,17),1] = 0.0 # add in element 21 if more than 5 profiles
  Lam_all[c(12,16,17,25),2] = 0.0 # add in element 21 if more than 5 profiles
  if (K == 6) {
    Lam_all[21,1] = 0.0
    Lam_all[21,2] = 0.0
    # If natural is profile 6:
    Lam_all[c(13,14,19,20,21), 6] = 0.0
    # If anthropogenic is profile 6:
    #Lam_all[c(order(Lam_all[1:(K-1),6]),15),6] = 0.0
    #Lam_all[order(Lam_all[1:(K-1),7]),7] = 0.0 # for an iron escape profile
  }
  if (K == 7) {
    Lam_all[c(8,11,12,17,19,21),1] = 0.0 # Baseline
    Lam_all[c(12,16,17,19,21,25),2] = 0.0 # Playa
    Lam_all[c(13,14,19,20,21,25),6] = 0.0 # Natural
    Lam_all[order(Lam_all[,7])[1:(K-1)],7] = 0.0 # Anthropogenic
  }
  
} else if (idf_rule == "something else") {
  
}
# ## and for the prior:
# for (k in 1:K ) {
#   alpha_matrix[k,which(Lam_all[,k] == 0.0)] = 5.0
#   beta_matrix[k,which(Lam_all[,k] == 0.0)] = 5.0*1e6
# }


############################################################################################################

# Make sure prior values are all assigned properly for the model in the list "dat_stan1":
dat_stan1 = list()

dat_stan1$L = dat_stan$L
dat_stan1$n = dat_stan$n
dat_stan1$alpha = t(alpha_matrix)
dat_stan1$beta = t(beta_matrix)
dat_stan1$flim = dat_stan$flim / 10000  ## base model analysis run on percent concentration (rather than ppm)
dat_stan1$n_obs = dat_stan$n_obs
dat_stan1$n_mis = dat_stan$n_mis
dat_stan1$n_bdl = dat_stan$n_bdl
dat_stan1$indx_obs = dat_stan$indx_obs
dat_stan1$indx_mis = dat_stan$indx_mis
dat_stan1$indx_bdl = dat_stan$indx_bdl
dat_stan1$y_obs = dat_stan$y_obs / 10000 ## base model analysis run on percent concentration (rather than ppm)
dat_stan1$dl_censor_vals = dat_stan$dl_censor_vals / 10000 ## base model analysis run on percent concentration (rather than ppm)
dat_stan1$K = K
dat_stan1$a = a
dat_stan1$b = b

dat_stan1$indx_non0Lam = which(Lam_all > 0.0) # vector indexing of the L by K Lambda matrix
dat_stan1$indx_0Lam = which(Lam_all == 0.0)
dat_stan1$n_non0Lam = length(dat_stan1$indx_non0Lam)
dat_stan1$n0Lam = length(dat_stan1$indx_0Lam)
all.equal(sort(c(dat_stan1$indx_non0Lam, dat_stan1$indx_0Lam)), 1:(dat_stan1$L*dat_stan1$K))



library("parallel")
n_chains =  4
(n_cores = min( detectCores(logical=FALSE), n_chains ))

# file1 = 'base_model.stan'
file1 = 'base_model_hard0.stan'

# Using the rstan library, run the model specified in "file1" on the data contained in dat_stan1.
# Use 8,000 iterations, discarding by default the first 2,000 as burn-in.
# 4 cores are used.
library("rstan")
fit = stan(
  file = file1, data = dat_stan1,
  chains = n_chains, cores = n_cores, iter = 8000, warmup = 2000,
  control = list(max_treedepth = 18,
                 adapt_delta = 0.999),
  pars = c("x", "lx", "lF", "lF0", "lF1", "llam_vec", "llam", "ly"), include=FALSE,
  seed = seed
)

# Save the file in an easily accessible location
save(fit, file=paste0('postsim/baseModel_hard0_idf_', idf_rule, '_K', K, '_inf', inflation_factor, "_seed", seed, '.rda'))
