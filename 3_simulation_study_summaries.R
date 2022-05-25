## This script should be run interactively and NOT sequentially or in batch mode.

rm(list=ls())
library("tidyverse")
library("readr")
library("lubridate")

files0 = list.files("postsim/")
length(files0)

# files = files0[grep("*se.*results_idfId2.*seed411", files0)] # Lam prior unimodal (alph >= 2); data cv ~ unif(0.2, 1.2); prior on cv mean .4 sd .15; 150 missing; adapt_delta = .95; Id2
# mesg = "Nmods_idId2_minalph2_seed411"

files = files0[grep("*eT.*results_K0.*seed418", files0)] # Lam prior unimodal (alph >= 2); data cv ~ unif(0.1, 0.8); prior on cv beta(1.2, 4.8); 150 missing; adapt_delta = .95; Id2
mesg = "sim_Tmods_idId2_minalph2_seed418"


(n_files = length(files))

sim_attr = lapply(files, function(x) {
  strsp = strsplit(x, "_")[[1]]
  data.frame(
    model = strsp[1],
    K0 = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("K0", strsp)]) ),
    Ke = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("Ke", strsp)]) ),
    src_used = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("present", strsp)]) ),
    infl = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("infl", strsp)]) ),
    sparsity = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("[[:digit:]]+sparse", strsp)]) ),
    seed = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("seed", strsp)]) )
  )
}) %>% do.call(rbind, .)
head(sim_attr)
tail(sim_attr)
dim(sim_attr)

sim_attrF = lapply(sim_attr[,1:ncol(sim_attr)], as.factor) %>% as.data.frame()
str(sim_attrF)

rhats_list = list()
f0_list = list()
Fmat_list = list()
Lam_prior_list = list()
Lam_truth_list = list()
mean_prof_list = list()
med_contrib_list = list()

FL_truth_list = list()
FL_hat_list = list()
lFL_resid_list = list()
lFL_nullresid_list = list()
FLpow_resid_list = list()
Y_list = list()
mseratios_list = list()
maeratios_list = list()

for (i in 1:n_files) {
  load(paste0("postsim/", files[i]))
  
  rhats_list[[i]] = rhats
  f0_list[[i]] = f_0
  Fmat_list[[i]] = F_mat
  Lam_prior_list[[i]] = Lam_prior
  Lam_truth_list[[i]] = Lam_truth
  mean_prof_list[[i]] = mean_profiles
  med_contrib_list[[i]] = as.matrix(median_contributions)
  
  FL_truth_list[[i]] = Lam_truth %*% rbind(f_0, F_mat) %>% t()  
  FL_hat_list[[i]] = exp(median_lFL)
  
  lFL_resid_list[[i]] = log(FL_truth_list[[i]]) - log(FL_hat_list[[i]])
  FLpow_resid_list[[i]] = (FL_truth_list[[i]])^0.2 - (FL_hat_list[[i]])^0.2
  
  Y_list[[i]] = t(Y)
    
  lFL_nullresid_list[[i]] = log(FL_truth_list[[i]]) - 
    log(tcrossprod(rep(1, nrow(FL_truth_list[[i]])), colMeans(Y_list[[i]], na.rm=TRUE)))
  
  mseratios_list[[i]] = mean(lFL_resid_list[[i]]^2) / mean(lFL_nullresid_list[[i]]^2)
  mseratios_list[[i]] = mean(abs(lFL_resid_list[[i]])) / mean(abs(lFL_nullresid_list[[i]]))
  
  cat(i, "of", n_files, "\r")
}

pairwisediff = function(var_vec, dfF=sim_attrF) {
  as.data.frame(cbind(var_vec, dfF)) %>% 
    pivot_wider(id_cols = c(K0, Ke, src_used, infl, sparsity, seed), names_from=model, values_from=var_vec) %>% 
    mutate(sparse = sparseT, base = baseT) %>%
    mutate(sparse_minus_base = sparse - base)
}

pairwiserat = function(var_vec, dfF=sim_attrF) {
  as.data.frame(cbind(var_vec, dfF)) %>% 
    pivot_wider(id_cols = c(K0, Ke, src_used, infl, sparsity, seed), names_from=model, values_from=var_vec) %>% 
    mutate(sparse = sparseT, base = baseT) %>%
    mutate(sparse_over_base = sparse / base)
}

paste_factor_names = function(src_used, sparsity, K0, Ke, infl, model=NULL) {
  if (is.null(model)) {
    paste0("src", src_used, " sp", sparsity, " K0", K0, " Ke", Ke, 
           " inf", sprintf("%1.1f", as.numeric(as.character(infl))))
  } else {
    paste0("src", src_used, " sp", sparsity, " K0", K0, " Ke", Ke, 
           " inf", sprintf("%1.1f", as.numeric(as.character(infl))), 
           " mod", toupper(substr(model, 1, 1)))
  }
}



### Convergence
rhats_list[[1]]
rhat_mean10worst = sapply(rhats_list, function(x) mean(x[order(x, decreasing=TRUE)[1:10]]))

## select one src_now and one sp_now
src_now = 1
src_now = 3
sp_now = 0
sp_now = 50

rhat_mean10worst_mod = lm(log(rhat_mean10worst) ~ K0*Ke*infl*model, 
                          data=as.data.frame(cbind(rhat_mean10worst, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now))
rhat_mean10worst_mod = lm(log(rhat_mean10worst) ~  Ke + infl + model + Ke:infl + Ke:model, 
                          data=as.data.frame(cbind(rhat_mean10worst, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now))
anova(rhat_mean10worst_mod)
summary(rhat_mean10worst_mod)
plot(rhat_mean10worst_mod)

ggplot(aes(x=rhat_mean10worst, 
           y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(rhat_mean10worst, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now) 
       ) + ylab("") + geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1))

ggsave(paste0("plots/rhats10worst_", mesg, "_src", src_now, "_sp", sp_now, ".pdf"), height=10, width=6)



library("xtable")

## percent of 20 highest Rhats exceeding 1.1
rhat_pct20worst = sapply(rhats_list, function(x) mean(x[order(x, decreasing=TRUE)[1:20]] > 1.1))


dat_now = as.data.frame(cbind(rhat_pct20worst, sim_attrF)) %>% 
  group_by(src_used, sparsity, model, K0, Ke, infl) %>% 
  summarize(Rhat=mean(rhat_pct20worst))
head(dat_now)  
tail(dat_now)
str(dat_now)

xt = xtabs(Rhat ~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
# xtc = xtabs(~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
(ft = ftable(round(xt, 3), col.vars = c(1,2,3), row.vars=c(4,5,6)))
xtableFtable(ft)



### Overall fit

mse_lLF = sapply(lFL_resid_list, function(x) mean(x^2))
mse_LFpow = sapply(FLpow_resid_list, function(x) mean(x^2))

mae_lLF = sapply(lFL_resid_list, function(x) mean(abs(x)))
mae_LFpow = sapply(FLpow_resid_list, function(x) mean(abs(x)))

mserat_lLF = sapply(mseratios_list, mean)
maerat_lLF = sapply(maeratios_list, mean)

## pickk one src_now and one sp_now value
src_now = 1
src_now = 3
sp_now = 0
sp_now = 50

mse_lLF_mod = lm(mse_lLF ~ K0*Ke*infl*model,
                 data=cbind(mse_lLF, sim_attrF) %>% filter(src_used==src_now, sparsity==sp_now))
mse_lLF_mod = lm(mse_lLF ~ Ke + infl + model,
                 data=cbind(mse_lLF, sim_attrF) %>% filter(src_used==src_now, sparsity==sp_now))

anova(mse_lLF_mod)
summary(mse_lLF_mod)
plot(mse_lLF_mod)


mserat_lLF_mod = lm(mserat_lLF ~ K0*Ke*infl*model,
                 data=cbind(mserat_lLF, sim_attrF) %>% filter(src_used==src_now, sparsity==sp_now))
mserat_lLF_mod = lm(mserat_lLF ~ K0 + Ke + infl + model + Ke:infl + K0:model + infl:model,
                 data=cbind(mserat_lLF, sim_attrF) %>% filter(src_used==src_now, sparsity==sp_now))

anova(mserat_lLF_mod)
summary(mserat_lLF_mod)
plot(mserat_lLF_mod)


msepwrat_lLF_mod = lm(sparse_over_base ~ K0*Ke*infl,
                 data=pairwiserat(mse_lLF) %>% filter(src_used==src_now, sparsity==sp_now))
msepwrat_lLF_mod = lm(sparse_over_base ~ K0 + Ke + infl + K0:Ke + Ke:infl,
                     data=pairwiserat(mse_lLF) %>% filter(src_used==src_now, sparsity==sp_now))
anova(msepwrat_lLF_mod)
summary(msepwrat_lLF_mod)
plot(msepwrat_lLF_mod)


ggplot(aes(x=log(mse_lLF), y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(mse_lLF, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now)
       ) + ylab("") + geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0))

ggsave(paste0("plots/mse_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=5)

ggplot(aes(x=mserat_lLF, y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(mserat_lLF, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now)
       ) + ylab("") + geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) 

ggsave(paste0("plots/mserat_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=5)

ggplot(aes(x=sparse_over_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwiserat(mse_lLF) %>% filter(src_used==src_now, sparsity==sp_now)
) + ylab("") + geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=1, lty=2) +
  ggtitle("MSE, pairwise ratio")

ggsave(paste0("plots/mse_pairwiserat_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=6, width=4)




library("xtable")

dat_pairs = pairwiserat(mse_lLF) %>% 
  group_by(sparsity, src_used, K0, Ke, infl) %>% 
  summarize(mse_lLF=mean(sparse_over_base), model="ratio")
dat_now = as.data.frame(cbind(mse_lLF,sim_attrF)) %>% 
  group_by(src_used, sparsity, model, K0, Ke, infl) %>% 
  summarize(mse_lLF=mean(mse_lLF)) %>%
  bind_rows(., dat_pairs[,c("src_used", "sparsity", "model", "K0", "Ke", "infl", "mse_lLF")])
dat_now$model = factor(dat_now$model, levels=c("baseT", "sparseT", "ratio"))
head(dat_now)  
tail(dat_now)
str(dat_now)



xt = xtabs(mse_lLF ~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
# xtc = xtabs(~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
(ft = ftable(round(xt, 3), col.vars = c(2,3), row.vars=c(1,4,5,6)))
xtableFtable(ft)



## Baseline profile
library("lattice")
i = 1
sim_attr[i,]
mean_prof_list[[i]]
Lam_truth_list[[i]]
levelplot(log(cbind(mean_prof_list[[i]][1,], Lam_truth_list[[i]][,1])+1e-8))


library("philentropy")

HellingerL1_dist = function(x, y) {
  stopifnot(  round(sum(x), 2) == round(sum(y), 2) & 
                round(sum(x), 2) == 1.0)
  sum( abs(sqrt(x) - sqrt(y)) )
}

pow = function(x, p) x^p - 1
pw = 0.4



### Unpaved dust (that is always present in simulation)
indx_K04 = which(sim_attr[,"K0"] == 4)
indx_unpaved = 4

JSD_unpaved_K04 = sapply(indx_K04, function(i) JSD(rbind(mean_prof_list[[i]][indx_unpaved,], Lam_truth_list[[i]][,"unpaved_alpha"]), unit="log", test.na=FALSE))
hist(JSD_unpaved_K04)

HellingerL1_unpaved_K04 = sapply(indx_K04, function(i) HellingerL1_dist(mean_prof_list[[i]][indx_unpaved,], Lam_truth_list[[i]][,"unpaved_alpha"]))
hist(HellingL1_unpaved_K04)

pow_unpaved_diff_K04 = sapply(indx_K04, function(i) sum(abs( pow(mean_prof_list[[i]][indx_unpaved,], pw) - pow(Lam_truth_list[[i]][,"unpaved_alpha"], pw) )))
hist(pow_unpaved_diff_K04)

## pick one src_now and one sp_now
src_now = 1
src_now = 3
sp_now = 0
sp_now = 50

ggplot(aes(x=log(JSD_unpaved_K04, base=10), y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(JSD_unpaved_K04, sim_attrF[indx_K04,])) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + 
  ggtitle("JSD Unpaved Dust profile")

ggsave(paste0("plots/JSD_unpaved_K04_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=4.5, width=5)


ggplot(aes(x=log(HellingerL1_unpaved_K04, base=10), y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(HellingerL1_unpaved_K04, sim_attrF[indx_K04,])) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + 
  ggtitle("Hellinger L1 \nUnpaved Dust profile")

ggsave(paste0("plots/HellingL1_unpaved_K04_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=4.5, width=5)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(log(JSD_unpaved_K04, base=2), sim_attrF[indx_K04,]) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + xlab("sparse minus base (log base 2)") +
  geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) + #+ xlim(-0.02, 0.02)
  ggtitle("JSD Unpaved Dust profile, pairwise difference\nK0=4")

ggsave(paste0("plots/JSD_unpaved_K04_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=6)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(log(HellingerL1_unpaved_K04, base=2), sim_attrF[indx_K04,]) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + xlab("sparse minus base (log base 2)") +
  geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) + #+ xlim(-0.02, 0.02)
  ggtitle("Hellinger L1 Unpaved Dust profile, pairwise difference\nK0=4")

ggsave(paste0("plots/JSD_unpaved_K04_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=6)



## Contributions
powF = 0.2 # 0.2 similar to log

mseFpw_unpaved_K04 = sapply(indx_K04, function(i) {
  unpaved_truth_indx = ifelse(sim_attr[i,"src_used"] == 1, 1, 3)
  mean( ( med_contrib_list[[i]][,indx_unpaved]^powF - Fmat_list[[i]][unpaved_truth_indx,]^powF )^2 )  
}  )

ggplot(aes(x=mseFpw_unpaved_K04, y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(mseFpw_unpaved_K04, sim_attrF[indx_K04,])) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) +
  ggtitle("MSE Unpaved Dust contributions\nK0=4")

ggsave(paste0("plots/MSEpowF_unpaved_K04_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=5, width=5)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(mseFpw_unpaved_K04, sim_attrF[indx_K04,])  %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("MSE (power) Unpaved Dust contributions\n pairwise difference, K0=4")

ggsave(paste0("plots/MSEpowF_unpaved_K04_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=6, width=7)


mseFpwK04_diff_mod = lm(sparse_minus_base ~ Ke*infl, 
                        data=pairwisediff(mseFpw_unpaved_K04, sim_attrF[indx_K04,]) %>% filter(src_used==src_now, sparsity==sp_now))

anova(mseFpwK04_diff_mod)
summary(mseFpwK04_diff_mod)

hist(med_contrib_list[[i]][,-c(indx_unpaved, 6:7)])
hist(med_contrib_list[[i]][,-c(indx_unpaved, 6:7)]^powF)

mseFpw_others = sapply(indx_K04, function(i) {
  if (sim_attr[i,"src_used"] == 1) {
    mean( (med_contrib_list[[i]][,2:(indx_unpaved-1)]^powF)^2 )
  } else if (sim_attr[i,"src_used"] == 3) {
    unpaved_truth_indx = 3
    mean( ( med_contrib_list[[i]][,2:(indx_unpaved-1)]^powF - t(Fmat_list[[i]][-unpaved_truth_indx,])^powF )^2 )
  }
})

ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(mseFpw_others, sim_attrF[indx_K04,])) + ylab("") + geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("MSE (power) other contributions\npairwise difference, K0=4")




## Unpaved dust profile (that is always present in simulation); comparing to "natural" profile
library("scales")

indx_K03 = which(sim_attr[,"K0"] == 3 & sim_attr[,"Ke"] > 0)
indx_natural = 4

JSD_unpaved_K03 = sapply(indx_K03, function(i) JSD(rbind(mean_prof_list[[i]][indx_natural,], Lam_truth_list[[i]][,"unpaved_alpha"]), unit="log", test.na=FALSE))
hist(JSD_unpaved_K03)

HellingerL1_unpaved_K03 = sapply(indx_K03, function(i) HellingerL1_dist(mean_prof_list[[i]][indx_natural,], Lam_truth_list[[i]][,"unpaved_alpha"]))
hist(HellingerL1_unpaved_K03)

pow_unpaved_diff_K03 = sapply(indx_K03, function(i) sum(abs( pow(mean_prof_list[[i]][indx_natural,], pw) - pow(Lam_truth_list[[i]][,"unpaved_alpha"], pw) )))
hist(pow_unpaved_diff_K03)


## pick one src_now and one sp_now
src_now = 1
src_now = 3
sp_now = 0
sp_now = 50

ggplot(aes(x=log(JSD_unpaved_K03, base=10), y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(JSD_unpaved_K03, sim_attrF[indx_K03,])) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + 
  ggtitle("JSD 'natural' profile\n(estimating Unpaved Dust)") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])
  
ggsave(paste0("plots/JSD_unpaved_K03_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=3.5, width=5)


ggplot(aes(x=log(HellingerL1_unpaved_K03, base=10), y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(HellingerL1_unpaved_K03, sim_attrF[indx_K03,])) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + 
  ggtitle("Hellinger L1 'natural' profile\n(estimating Unpaved Dust)") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])

ggsave(paste0("plots/HellingL1_unpaved_K03_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=3.5, width=5)


ggplot(aes(x=pow_unpaved_diff_K03, y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke),
       data=as.data.frame(cbind(pow_unpaved_diff_K03, sim_attrF[indx_K03,])) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + 
  scale_fill_manual(values=hue_pal()(3)[2:3])


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(log(JSD_unpaved_K03, base=2), sim_attrF[indx_K03,]) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + xlab("sparse minus base (log base 2)") + 
  geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("JSD 'natural' profile (estimating Unpaved Dust)\npairwise difference, K0=3") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])

ggsave(paste0("plots/JSD_unpaved_K03_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=6, width=7)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(log(HellingerL1_unpaved_K03, base=2), sim_attrF[indx_K03,]) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + xlab("sparse minus base (log base 2)") + 
  geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("Hellinger L1 'natural' profile (estimating Unpaved Dust)\npairwise difference, K0=3") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])

ggsave(paste0("plots/HellingL1_unpaved_K03_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=6, width=7)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(pow_unpaved_diff_K03, sim_attrF[indx_K03,]) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("Power loss 'natural' profile (estimating Unpaved Dust)\npairwise difference, K0=3") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])



## contributions
powF = 0.2 # 0.2 similar to log

mseFpw_unpaved_K03 = sapply(indx_K03, function(i) {
  unpaved_truth_indx = ifelse(sim_attr[i,"src_used"] == 1, 1, 3)
  mean( ( med_contrib_list[[i]][,indx_natural]^powF - Fmat_list[[i]][unpaved_truth_indx,]^powF )^2 )  
}  )

ggplot(aes(x=mseFpw_unpaved_K03, y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(mseFpw_unpaved_K03, sim_attrF[indx_K03,])) %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) +
  ggtitle("MSE 'natural' contributions\n(estimating Unpaved Dust), K0=3") + 
  scale_fill_manual(values=hue_pal()(3)[2:3])

ggsave(paste0("plots/MSEpowF_unpaved_K03_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=3.75, width=5)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(mseFpw_unpaved_K03, sim_attrF[indx_K03,])  %>% filter(src_used==src_now, sparsity==sp_now)) + 
  ylab("") + geom_boxplot(show.legend=FALSE) + geom_vline(xintercept=0, lty=2) +
  ggtitle("MSE (power) 'natural' contributions (estimating Unpaved Dust)\npairwise difference, K0=3",) + 
  scale_fill_manual(values=hue_pal()(3)[2:3])

ggsave(paste0("plots/MSEpowF_unpaved_K03_pairwisediff_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=6, width=7)


mseFpwK03_diff_mod = lm(sparse_minus_base ~ Ke*infl, 
                        data=pairwisediff(mseFpw_unpaved_K03, sim_attrF[indx_K03,]) %>% filter(src_used==src_now, sparsity==sp_now))

anova(mseFpwK03_diff_mod)
summary(mseFpwK03_diff_mod)


### discrimination between present and absent sources
Fpw_discrim = sapply(1:nrow(sim_attr), function(i) 
  # for(i in 1:nrow(sim_attr))
  {
  if (sim_attr[i,"src_used"] == 1) {
    if (sim_attr[i,"K0"] == 3) {
      if (sim_attr[i, "Ke"] > 0) {
          kindx_active = 4 # natural
          kindx_active_truth = 1
          cancompute = TRUE
      } else {
        cancompute = FALSE
      }
    } else if (sim_attr[i,"K0"] == 4) {
        kindx_active = 4 # unpaved
        kindx_active_truth = 1
        cancompute = TRUE
    }
  } else if (sim_attr[i,"src_used"] == 3) {
    if (sim_attr[i,"K0"] == 3) {
      if (sim_attr[i, "Ke"] > 1) {
          kindx_active = 2:4 # two sources and natural
          kindx_active_truth = 1:3
          cancompute = TRUE
      } else {
        cancompute = FALSE
      }
    } else if (sim_attr[i,"K0"] == 4) {
      if (sim_attr[i, "Ke"] > 1) {
          kindx_active = 2:4 # two sources and natural
          kindx_active_truth = 1:3
          cancompute = TRUE
      } else {
        cancompute = FALSE
      }
    }
  }

  if (cancompute) {
    if (sim_attr[i,"sparsity"] == 0) {
      med_contrib_active = med_contrib_list[[i]][,kindx_active]
    } else {
      nonzero_active_truth = which( t(Fmat_list[[i]][kindx_active_truth,]) > 0 )
      med_contrib_active = med_contrib_list[[i]][,kindx_active][nonzero_active_truth]
    }
    med_contrib_inactive = med_contrib_list[[i]][,-c(1,kindx_active)]
  } 
  
  if (cancompute) {
    out = (quantile(med_contrib_active, 0.75) / quantile(med_contrib_inactive, 0.75))^0.5
  } else {
    out = NA
  }
  out
  # cat(i, "of", nrow(sim_attr), "\r")
}
)


## pick one src_now and one sp_now
src_now = 1
src_now = 3
sp_now = 0
sp_now = 50

ggplot(aes(x=Fpw_discrim, y=paste_factor_names(src_used, sparsity, K0, Ke, infl, model), fill=Ke, alpha=1.0*grepl("base", model)),
       data=as.data.frame(cbind(Fpw_discrim, sim_attrF)) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + geom_vline(xintercept=1, lty=2) + 
  ggtitle("Fpow ratio active/inactive")

ggsave(paste0("plots/Fpow_discrim_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=6)


ggplot(aes(x=sparse_minus_base, y=paste_factor_names(src_used, sparsity, K0, Ke, infl), fill=Ke),
       data=pairwisediff(Fpw_discrim, sim_attrF) %>% filter(src_used==src_now, sparsity==sp_now)) + ylab("") + 
  geom_boxplot(show.legend=FALSE) + scale_alpha(range=c(0.3, 1.0)) + geom_vline(xintercept=0, lty=2) + 
  ggtitle("Fpow ratio active/inactive\nPairwise difference")

ggsave(paste0("plots/Fpow_discrim_pairwise_src", src_now, "_sp", sp_now, "_", mesg, ".pdf"), height=8, width=6)



## Tables
library("xtable")

dat_pairs = pairwiserat(HellingerL1_unpaved_K03, sim_attrF[indx_K03,]) %>% 
  bind_rows(., pairwiserat(HellingerL1_unpaved_K04, sim_attrF[indx_K04,])) %>%
  group_by(sparsity, src_used, K0, Ke, infl) %>% 
  summarize(HellingerL1_unpaved=mean(sparse_over_base), model="ratio")
dat_now = as.data.frame(cbind(HellingerL1_unpaved=HellingerL1_unpaved_K03, sim_attrF[indx_K03,])) %>% 
  bind_rows(., cbind(HellingerL1_unpaved=HellingerL1_unpaved_K04, sim_attrF[indx_K04,])) %>%
  group_by(src_used, sparsity, model, K0, Ke, infl) %>% 
  summarize(HellingerL1_unpaved=mean(HellingerL1_unpaved)*10) %>%
  bind_rows(., dat_pairs[,c("src_used", "sparsity", "model", "K0", "Ke", "infl", "HellingerL1_unpaved")])
dat_now$model = factor(dat_now$model, levels=c("baseT", "sparseT", "ratio"))
head(dat_now)  
tail(dat_now)
# str(dat_now)

xt = xtabs(HellingerL1_unpaved ~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
# xtc = xtabs(~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
(ft = ftable(round(xt, 2), col.vars = c(2,3), row.vars=c(1,4,5,6)))
xtableFtable(ft)



dat_pairs = pairwiserat(mseFpw_unpaved_K03, sim_attrF[indx_K03,]) %>% 
  bind_rows(., pairwiserat(mseFpw_unpaved_K04, sim_attrF[indx_K04,])) %>%
  group_by(sparsity, src_used, K0, Ke, infl) %>% 
  summarize(mseFpw_unpaved=mean(sparse_over_base), model="ratio")
dat_now = as.data.frame(cbind(mseFpw_unpaved=mseFpw_unpaved_K03, sim_attrF[indx_K03,])) %>% 
  bind_rows(., cbind(mseFpw_unpaved=mseFpw_unpaved_K04, sim_attrF[indx_K04,])) %>%
  group_by(src_used, sparsity, model, K0, Ke, infl) %>% 
  summarize(mseFpw_unpaved=mean(mseFpw_unpaved)*1) %>%
  bind_rows(., dat_pairs[,c("src_used", "sparsity", "model", "K0", "Ke", "infl", "mseFpw_unpaved")])
dat_now$model = factor(dat_now$model, levels=c("baseT", "sparseT", "ratio"))
head(dat_now)  
tail(dat_now)
# str(dat_now)

xt = xtabs(mseFpw_unpaved ~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
# xtc = xtabs(~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
(ft = ftable(round(xt, 2), col.vars = c(2,3), row.vars=c(1,4,5,6)))
xtableFtable(ft)



dat_pairs = pairwisediff(Fpw_discrim) %>% 
  group_by(sparsity, src_used, K0, Ke, infl) %>% 
  summarize(Fpw_discrim=mean(sparse_minus_base), model="diff")
dat_now = as.data.frame(cbind(Fpw_discrim, sim_attrF)) %>% 
  group_by(src_used, sparsity, model, K0, Ke, infl) %>% 
  summarize(Fpw_discrim=mean(Fpw_discrim)) %>%
  bind_rows(., dat_pairs[,c("src_used", "sparsity", "model", "K0", "Ke", "infl", "Fpw_discrim")])
dat_now$model = factor(dat_now$model, levels=c("baseT", "sparseT", "diff"))
head(dat_now)  
tail(dat_now)
str(dat_now)

xt = xtabs(Fpw_discrim ~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
# xtc = xtabs(~ src_used + sparsity + model + K0 + Ke + infl, data=dat_now)
(ft = ftable(round(xt, 3), col.vars = c(2,3), row.vars=c(1,4,5,6)))
xtableFtable(ft)











### profile plots
## K04
indx_unpaved = 4

pow = function(x, p) x^p - 1
pw = 0.4
ticks = c(0, 0.01, 0.1, 0.25, 0.5, 1.0)
elems_list <- rownames(Lam_truth)

pdf(file=paste0("plots/profiles_", mesg, "_k0", 4, ".pdf"), width=8, height=5)
  for (src in unique(sim_attr$src_used)) {
    for (spc in unique(sim_attr$sparsity)) {
      for (inf in unique(sim_attr$infl)) {
        
      indx_now = which(sim_attr$K0 == 4 & sim_attr$src_used == src & sim_attr$infl == inf & sim_attr$sparsity == spc)
      
      plot(pow(Lam_truth_list[[indx_now[[1]]]][,"unpaved_alpha"], pw) + 1, axes=F, ylab="", xlab="", ylim=c(0,1), pch=19)
      axis(side=1, at=1:length(elems_list), labels = elems_list)
      axis(side=2, at=1+pow(ticks,pw), labels=ticks)
      
      for (i in indx_now) {
        points(pow(mean_prof_list[[i]][indx_unpaved,], pw)+1, 
               col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
               pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4))
      }
      points(pow(Lam_truth_list[[indx_now[1]]][,"unpaved_alpha"], pw) + 1, pch=19, cex=1.5)
      title(main=paste0("Unpaved Dust profile\nsrc=", src, " spc=", spc, "\nK0=4 inf=", inf))
      legend("topright", pch=c(3,4,19), col=c("blue", "red", "black"), legend=c("sparse", "base", "truth"), bty="n")
    } 
  }
}
dev.off()


## K03

indx_natural = 4

pdf(file=paste0("plots/profiles_", mesg, "_k0", 3, ".pdf"), width=11, height=7)
  for (src in unique(sim_attr$src_used)) {
    for (spc in unique(sim_attr$sparsity)) {
      for (ke in 1:2) {
        for (inf in unique(sim_attr$infl)) {
        
        indx_now = which(sim_attr$K0 == 3 & sim_attr$Ke == ke & sim_attr$src_used == src & sim_attr$infl == inf & sim_attr$sparsity == spc)
        
        plot(pow(Lam_truth_list[[indx_now[1]]][,"unpaved_alpha"], pw) + 1, axes=F, ylab="", xlab="", ylim=c(0,1), pch=19)
        axis(side=1, at=1:length(elems_list), labels = elems_list)      
        axis(side=2, at=1+pow(ticks,pw), labels=ticks)
        
        for (i in indx_now) {
          points(pow(mean_prof_list[[i]][indx_natural,], pw)+1, 
                 col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
                 pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4))
        }
        points(pow(Lam_truth_list[[indx_now[1]]][,"unpaved_alpha"], pw) + 1, pch=19, cex=1.5)
        points(pow(Lam_prior_list[[indx_now[1]]][,indx_natural], pw) + 1, pch=1, cex=1.5)
        points(pow(Lam_truth_list[[indx_now[1]]][,"baseline_alpha"], pw) + 1, pch=17)
        title(main=paste0("Unpaved Dust profile (estimated with 'natural')\nsrc=", src, " spc=", spc, "\nK0=3 inf=", inf, " Ke=", ke))
        legend("topright", pch=c(3,4,19,1,17), col=c("blue", "red", "black", "black", "black"), legend=c("sparse", "base", "unpaved truth", "natural prior", "baseline truth"), bty="n")
      } 
    }
  }  
}

dev.off()




### contrib plots
## K04
powF = 0.2 # 0.2 similar to log

indx_unpaved = 4 # can change this from 4 to look at other model contribs, just ignore "truth"
ticksF = 10.0^(0:5)
nloc = 30

pdf(file=paste0("plots/contrib_", mesg, "_k0", 4, ".pdf"), width=11, height=7)
par(mar=c(5,4,5,1)+0.1)
for (src in unique(sim_attr$src_used)) {
  for (spc in unique(sim_attr$sparsity)) {
    for (inf in unique(sim_attr$infl)) {
      for(seed in unique(sim_attr$seed)[1:6]) {

        indx_now = which(sim_attr$K0 == 4 & sim_attr$src_used == src & sim_attr$infl == inf & sim_attr$sparsity == spc & sim_attr$seed == seed)
        
        if (length(indx_now)>0) {
          unpaved_truth_indx = ifelse(sim_attr[indx_now[1],"src_used"] == 1, 1, 3)
          rhatsF_now = sapply(indx_now, function(ii) mean(rhats_list[[ii]][grep(paste0("F\\[", indx_unpaved, ","), names(rhats_list[[ii]]))][1:nloc]))
          
          plot(Fmat_list[[indx_now[1]]][unpaved_truth_indx, 1:nloc]^powF, axes=F, ylab="ppm", xlab="", ylim=c(0,15), pch=19)
          axis(side=1)      
          axis(side=2, at=ticksF^powF, labels=ticksF)
          
          for (i in indx_now) {
            points(med_contrib_list[[i]][1:nloc, indx_unpaved]^powF, 
                   col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
                   # pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4)
                   pch=as.character(sim_attr[i, "Ke"])
                   )
          }
          points(Fmat_list[[indx_now[1]]][unpaved_truth_indx, 1:nloc]^powF, pch=19, cex=1.0)
          title(main=paste0("Unpaved Dust contributions\nsrc=", src, " spc=", spc, "\nK0=4 inf=", inf, 
                            "\nseed=", seed))
          legend("topright", pch=c(3,4), col=c("blue", "red"), legend=c("sparse", "base"), bty="n") 
          legend("topleft", ncol=4, title="Rhat", #col=rep(c("red", "blue"), times=3),
                 legend=c("base", "sparse", round(tapply(rhatsF_now, sim_attr[indx_now, c("model", "Ke")], mean), 2)), bty="n")
        }
      }
    } 
  }
}

dev.off()


## K03 natural contrib vs unpaved contrib
indx_natural = 4 # can change this from 4 to look at other model contribs, just ignore "truth"

pdf(file=paste0("plots/contrib_", mesg, "_k0", 3, ".pdf"), width=11, height=7)
par(mar=c(5,4,5,1)+0.1)
for (src in unique(sim_attr$src_used)) {
  for (spc in unique(sim_attr$sparsity)) {
    # for (ke in 1:2) {
      for (inf in unique(sim_attr$infl)) {
        for(seed in unique(sim_attr$seed)[1:2]) {
          
          indx_now = which(sim_attr$K0 == 3 & sim_attr$src_used == src & sim_attr$infl == inf & sim_attr$sparsity == spc & sim_attr$seed == seed & sim_attr$Ke > 0 )
          unpaved_truth_indx = ifelse(sim_attr[indx_now[1],"src_used"] == 1, 1, 3)
          rhatsF_now = sapply(indx_now, function(ii) mean(rhats_list[[ii]][grep(paste0("F\\[", indx_natural, ","), names(rhats_list[[ii]]))][1:nloc]))
          
          plot(Fmat_list[[indx_now[1]]][unpaved_truth_indx, 1:nloc]^powF, axes=F, ylab="ppm", xlab="", ylim=c(0,15), pch=19)
          axis(side=1)
          axis(side=2, at=ticksF^powF, labels=ticksF)
          
          for (i in indx_now) {
            points(med_contrib_list[[i]][1:nloc, indx_natural]^powF, 
                   col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
                   # pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4)
                   pch=as.character(sim_attr[i, "Ke"])
                   )
          }
          points(Fmat_list[[indx_now[1]]][unpaved_truth_indx, 1:nloc]^powF, pch=19, cex=1.0)
          title(main=paste0("Unpaved Dust contributions (estimated w Natural)\nsrc=", src, " spc=", spc, 
                            "\nK0=3 inf=", inf,  "\nseed=", seed))
          legend("topright", pch=c(3,4),
                 col=c("blue", "red"), legend=c("sparse", "base"), bty="n")
          legend("topleft", ncol=3, title="Rhat", #col=rep(c("red", "blue"), times=3),
                 legend=c("base", "sparse", round(tapply(rhatsF_now, sim_attr[indx_now, c("model", "Ke")], mean), 2)), bty="n")
        }
      }
    # }
  }
}
dev.off()

