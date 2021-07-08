# Load in the raw rhizoplaca lichen data
load(file = 'data/dat_rhizICP.rda')

# Run the file below to access prior information for each profile
source('helperScripts/all_profiles_nitrogen_modified_specificInflation.R')

# Assign prior
L = dat_stan$L

## prior on error coefficient of variation
mean_cv = rep(0.5, L)
sd_cv = rep(0.2, L)

a = mean_cv^2 / sd_cv^2
b = mean_cv / sd_cv^2

# Use baseline, playa, brake, exhaust, and unpaved profiles.
# Impose zeroes in specific places in order to ensure identifiability conditions are satisfied.
src_use = c(1,2,4,5,6)
alpha_matrix = rbind(alpha_matrix[src_use,])
beta_matrix = rbind(beta_matrix[src_use,])
el_use = 1:ncol(alpha_matrix)
(K = nrow(alpha_matrix))
(L = length(el_use))
x = t(alpha_matrix[, el_use] / beta_matrix[, el_use])
(Lam_all = x / rep(colSums(x), each=nrow(x)))
for (k in 3:5 ) {
  ord = order(Lam_all[,k])
  Lam_all[ord[1:(K-1)],k] = 1e-9
}
Lam_all[c(8,11,12,17),1] = 1e-9 # add in element 21 if more than 5 profiles
Lam_all[c(12,16,17,25),2] = 1e-9 # add in element 21 if more than 5 profiles
# If natural is profile 6:
#Lam_all[c(13,14,19,20,21),6] = 1e-9
# If anthropogenic is profile 6:
#Lam_all[c(order(Lam_all[1:(K-1),6]),15),6] = 1e-9
#Lam_all[order(Lam_all[1:(K-1),7]),7] = 1e-9 # for an iron escape profile

# If 7 profiles: ##################################
#Lam_all[c(8,11,12,17,19,21),1] = 1e-9 # Baseline
#Lam_all[c(12,16,17,19,21,25),2] = 1e-9 # Playa
# Lam_all[c(13,14,19,20,21,25),6] = 1e-9 # Natural
# Lam_all[order(Lam_all[1:(K-1),7]),7] = 1e-9 # Anthropogenic

Lam_all
## and for the prior:
for (k in 1:K ) {
  alpha_matrix[k,which(Lam_all[,k] == 1e-9)] = 5.0
  beta_matrix[k,which(Lam_all[,k] == 1e-9)] = 5.0*1e6
}
#############################################################################################################

# Make sure prior values are all assigned properly for the model:
dat_stan$alpha = t(alpha_matrix)
dat_stan$beta = t(beta_matrix)

dat_stan$a = a
dat_stan$b = b
dat_stan$K = K

library("parallel")
n_chains =  4
(n_cores = min( detectCores(logical=FALSE), n_chains ))

file1 = 'sparse_model_margF.stan'

# Using the rstan library, run the model specified in "file1" on the data contained in dat_stan.
# Use 10,000 iterations, discarding by default the first 5,000 as burn-in.
# 4 cores are used.
library("rstan")
fit = stan(
  file = file1, data = dat_stan,
  chains = n_chains, cores = n_cores, iter = 10000,
  control = list(max_treedepth = 18,
                 adapt_delta = 0.9999),
  pars = c("x", "lx", "lF", "lF0", "llam"), include=FALSE
)

# Save the file in an easily accessible location
save(fit, file=paste0('postsim/sparseModel_', K, 'prof_identified', '.rda'))
