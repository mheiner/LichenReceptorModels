rm(list=ls())

library("loo")
library("tidyverse")


files0 = list.files("postsim/")
length(files0)

# files = files0[grep("fit.*null.*Model.*_idf_Id2.*seed427", files0)] # production runs
# mesg = "null_mods_idId2_minalph5_seed427"

# files = files0[grep("fit.*eModel.*_idf_Id2.*seed427", files0)] # production runs with normal likelihood
# mesg = "Nmods_idId2_minalph5_seed427"

files = files0[grep("fit.*eTModel.*_idf_Id2.*seed427", files0)] # production runs with t likelihood
mesg = "Tmods_idId2_minalph5_seed427"

(n_files = length(files))
files

for (i in 1:n_files) {
  load(paste0("postsim/", files[i]))
  if (any(grepl("llik", names(fit)))) {
    llik = extract_log_lik(fit, parameter_name = "llik", merge_chains=FALSE)
    
    r_eff = relative_eff(exp(llik), cores=2)
    summary(r_eff)
    loo_fit = loo(llik, r_eff=r_eff, cores=4)
    cat(files[i],"\n")
    print(loo_fit)
  }
}
