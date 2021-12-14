rm(list=ls())

# source("all_profiles.R")
# source("all_profiles_nitrogen_modified.R")
inflation_factor = 3.0
source("helperScripts/all_profiles_nitrogen_modified_specificInflation.R")

ls()
dim(alpha_matrix)
dim(beta_matrix)

(src_names = unname(sapply(rownames(alpha_matrix), function(x) strsplit(x, "_")[[1]][1])))
(el_names = colnames(alpha_matrix))

src_use = 1:nrow(alpha_matrix) # if using all profiles
src_use = c(1:5) # or select a subset
# (K = length(src_use))

el_use = 1:ncol(alpha_matrix)
el_use = c(1, 3, 4, 7, 15, 18, 22, 23)
(L = length(el_use))

sim_Lam_prior = function(alpha_matrix, beta_matrix) {
  X = matrix(rgamma(prod(dim(alpha_matrix)), alpha_matrix, rate=beta_matrix), nrow=nrow(alpha_matrix))
  sums = rowSums(X)
  Out = X / sums
  return(Out)
}

n = 1e4
Lsim = replicate(n, sim_Lam_prior(alpha_matrix, beta_matrix))
dim(Lsim)
dimnames(Lsim) = list(src_names, el_names, NULL)

Lsim_mean = apply(Lsim, c(1,2), mean)
Gam_mean = alpha_matrix / beta_matrix
library("lattice")
levelplot(Lsim_mean[src_use,])
levelplot(Gam_mean / rowSums(Gam_mean))

hist(Lsim[2,1,])
hist(Lsim[4,14,])

src_names

pdf(file=paste0("infF", inflation_factor, ".pdf"), height=20, width=6)
for (k in src_use) {
  # k = 4
  par(mfrow=c(L,1))
  for (j in el_use) {
    hist(Lsim[k,j,], xlim=c(0,1), main=paste0(src_names[k], ": ", el_names[j]))
  }
}
dev.off()



par(mfrow=c(7,1))
k = 6
hist(Lsim[k,1,], xlim=c(0,1))
hist(Lsim[k,2,], xlim=c(0,1))
hist(Lsim[k,4,], xlim=c(0,1))
hist(Lsim[k,7,], xlim=c(0,1))
hist(Lsim[k,15,], xlim=c(0,1))
hist(Lsim[k,18,], xlim=c(0,1))
hist(Lsim[k,22,], xlim=c(0,1))

colnames(alpha_matrix)


# nitrogen
hist(Lsim[1,4,])
hist(Lsim[2,4,])
hist(Lsim[3,4,])
hist(Lsim[4,4,])

hist(Lsim[6,4,]) # unpaved, Nitrogen
summary(Lsim[6,4,])
dim(Lsim)

alpha_matrix / beta_matrix
sqrt(alpha_matrix / beta_matrix^2)[6,]
