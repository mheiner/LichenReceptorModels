ls()
dim(alpha_matrix)
dim(beta_matrix)

(src_names = unname(sapply(rownames(alpha_matrix), function(x) strsplit(x, "_")[[1]][1])))
(el_names = colnames(alpha_matrix))

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

Lsim_med = apply(Lsim, c(1,2), median)
Lsim_mean = apply(Lsim, c(1,2), mean)
Lsim_sd = apply(Lsim, c(1,2), sd)

Lsim_sd / Lsim_mean
