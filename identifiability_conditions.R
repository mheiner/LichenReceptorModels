## sensitivity settings
rm(list=ls())

K = 5
K = 6
K = 7
inflation_factor = 1.0 # replaces the assignment at beginning of "all_prof...R"


# Load in the raw rhizoplaca lichen data
load(file = 'data/dat_rhizICP.rda')

# Run the file below to access prior information for each profile
source('helperScripts/all_profiles_nitrogen_modified_specificInflation.R')
L = dat_stan$L


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
(K = nrow(alpha_matrix))
(L = length(el_use))

source('helperScripts/simProfiles_prior.R')
dim(Lsim)
dimnames(Lsim)
Lsim_med = t(Lsim_med)
Lsim_sd = t(Lsim_sd)

x = t(alpha_matrix[, el_use] / beta_matrix[, el_use])

Lam = x / rep(colSums(x), each=nrow(x))
colnames(Lam)
colnames(Lsim_med)
colnames(Lam) = colnames(Lsim_med)

Lam_orig = Lam

## original identifiability conditions
for (k in 3:5) {
  ord = order(Lam_orig[,k])
  Lam_orig[ord[1:(K-1)],k] = 1e-9
}
Lam_orig[c(8,11,12,17),1] = 1e-9 # add in element 21 if more than 5 profiles
Lam_orig[c(12,16,17,25),2] = 1e-9 # add in element 21 if more than 5 profiles
if (K == 6) {
  Lam_orig[21,1] = 1e-9
  Lam_orig[21,2] = 1e-9
  # If natural is profile 6:
  Lam_orig[c(13,14,19,20,21), 6] = 1e-9
  # If anthropogenic is profile 6:
  #Lam_orig[c(order(Lam_orig[1:(K-1),6]),15),6] = 1e-9
  #Lam_orig[order(Lam_orig[1:(K-1),7]),7] = 1e-9 # for an iron escape profile
}
if (K == 7) {
  Lam_orig[c(8,11,12,17,19,21),1] = 1e-9 # Baseline
  Lam_orig[c(12,16,17,19,21,25),2] = 1e-9 # Playa
  Lam_orig[c(13,14,19,20,21,25),6] = 1e-9 # Natural
  Lam_orig[order(Lam_orig[,7])[1:(K-1)],7] = 1e-9 # Anthropogenic
}

pow = function(x, p) x^p
pw = 0.4

library("lattice")
# pdf(file="Lam_prior_means.pdf", width=8, height=4)
levelplot(pow(Lam[,K:1], pw), main="Lam")
# dev.off()

# pdf(file="Lam_sim_medians.pdf", width=8, height=4)
levelplot(pow(Lsim_med[,K:1], pw), main="Lam sim median")
# dev.off()

# pdf(file="Lam_prior_sds.pdf", width=8, height=4)
levelplot(pow(Lsim_sd[,K:1], pw), main="Lam sim st. deviation")
# dev.off()

Z = Lam_orig # from model spec

dim(Z)
(K = ncol(Z))
Z[Z==1e-9] = 0
Z

# pdf(file="Lam_id01.pdf", width=8, height=4)
levelplot(1*(Z[,K:1] > 0.0), main="Lam identified 1-0", col.regions=gray((0:15)/15))
# dev.off()