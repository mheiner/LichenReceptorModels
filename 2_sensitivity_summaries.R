## This script should be run interactively and NOT sequentially or in batch mode

rm(list=ls())
library("tidyverse")
library("readr")
library("lubridate")

files0 = list.files("postsim/")
length(files0)



# files = files0[grep("results.*TModel.*_idf_Id2.*seed418", files0)] # Lam prior unimodal (alph > 2); beta(1.2, 4.8) prior on "cv"; adapt_delta = .95; K0 4, 5
# mesg = "Tmods_idId2_minalph2_seed418"

# files = files0[grep("results.*TModel.*_idf_Id2.*seed419", files0)] # Lam prior unimodal (alph > 5); beta(1.2, 4.8) prior on "cv"; adapt_delta = .95; K0 4, 5
# mesg = "Tmods_idId2_minalph5_seed419"

# files = files0[grep("results.*TModel.*_idf_original.*seed419", files0)] # Lam prior unimodal (alph > 5); beta(1.2, 4.8) prior on "cv"; adapt_delta = .95; K0 5 only
# mesg = "Tmods_idOrig_minalph5_seed419"


# files = files0[grep("results.*null.*Model.*_idf_Id2.*seed427", files0)] # production runs
# mesg = "Nmods_idId2_minalph5_seed427"

# files = files0[grep("results.*eModel.*_idf_Id2.*seed427", files0)] # production runs
# mesg = "Nmods_idId2_minalph5_seed427"

files = files0[grep("results.*eTModel.*_idf_Id2.*seed427", files0)] # production runs
mesg = "Tmods_idId2_minalph5_seed427"



(n_files = length(files))

sim_attr = lapply(files, function(x) {
  strsp = strsplit(x, "_")[[1]]
  data.frame(
    model = gsub("Model", "", strsp[2]),
    K0 = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("K0", strsp)]) ),
    Ke = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("Ke", strsp)]) ),
    infl = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("inf", strsp)]) ),
    seed = as.numeric( gsub("[[:alpha:]]", "", strsp[grep("seed", strsp)]) )
  )
}) %>% do.call(rbind, .)
head(sim_attr)
tail(sim_attr)
dim(sim_attr)

sim_attrF = lapply(sim_attr[,1:ncol(sim_attr)], as.factor) %>% as.data.frame()
str(sim_attrF)
table(sim_attrF[,1:4])

rhats_list = list()
Lam_prior_list = list()
mean_prof_list = list()
med_contrib_list = list()

FL_hat_list = list()
resid_list = list()
nullresid_list = list()
pow_resid_list = list()
Y_list = list()
mseratios_list = list()
maeratios_list = list()

for (i in 1:n_files) {
  load(paste0("postsim/", files[i]))
  
  rhats_list[[i]] = rhats[grep("cv|^F|lam|nu", names(rhats))]
  Lam_prior_list[[i]] = dat_stan$Lam_prior / rep(colSums(dat_stan$Lam_prior), each=nrow(dat_stan$Lam_prior))
  
  if (!grepl("null", sim_attr[i,"model"])) {
    median_contributions <- summ$summary[ grep("^F", rownames(summ$summary)) , "50%"]
    dim(median_contributions) <- c(dat_stan$n, dat_stan$K)
    median_contributions <- as.data.frame(median_contributions)
    
    lFL = summ$summary[grep("lLF", rownames(summ$summary)), "50%"]# arranged with sample index changing first
    dim(lFL) = c(dat_stan$n, dat_stan$L)
  }
  
  mean_profiles <- summ$summary[ grep("lam", rownames(summ$summary)) , "mean"]
  dim(mean_profiles) <- c(dat_stan$K, dat_stan$L)
  
  stopifnot( max(abs(rowSums(mean_profiles) - 1.0)) < 1e-14 )
  
  cvs <- summ$summary[ grep("cv", rownames(summ$summary)) , "mean"]

  Y0 = matrix(NA, nrow=dat_stan$L, ncol=dat_stan$n)
  Y0[dat_stan$indx_obs] = dat_stan$y_obs
  Y = t(Y0)
  rm(Y0)
  
  Y_list[[i]] = Y
  
  mean_prof_list[[i]] = mean_profiles
  
  if (!grepl("null", sim_attr[i,"model"])) {
    FL_hat_list[[i]] = exp(lFL)
 
    med_contrib_list[[i]] = as.matrix(median_contributions)
  
    resid_list[[i]] = log(Y_list[[i]]) - log(FL_hat_list[[i]])
    pow_resid_list[[i]] = (Y_list[[i]])^0.2 - (FL_hat_list[[i]])^0.2
  
    nullresid_list[[i]] = log(Y_list[[i]]) - 
      log(tcrossprod(rep(1, nrow(Y_list[[i]])), colMeans(Y_list[[i]], na.rm=TRUE)))
  
    mseratios_list[[i]] = (apply(resid_list[[i]]^2, 2, mean, na.rm=TRUE)) / (apply(nullresid_list[[i]]^2, 2, mean, na.rm=TRUE))
    maeratios_list[[i]] = (apply(abs(resid_list[[i]]), 2, mean, na.rm=TRUE)) / (apply(abs(nullresid_list[[i]]), 2, mean, na.rm=TRUE))
  }
  
  cat(i, "of", n_files, "\r")
}

pairwisediff = function(var_vec, dfF=sim_attrF) {
  as.data.frame(cbind(var_vec, dfF)) %>% 
    pivot_wider(id_cols = c(K0, Ke, infl, seed), names_from=model, values_from=var_vec) %>% 
    mutate(sparse = sparseT, base = baseT) %>%
    mutate(sparse_minus_base = sparse - base)
}


## Convergence
head(rhats_list[[1]])

rhat_mean10worst = sapply(rhats_list, function(x) mean(x[order(x, decreasing=TRUE)[1:10]]))
hist(rhat_mean10worst)

ggplot(aes(x=rhat_mean10worst, y=paste0("K0", K0, "_Ke", Ke, "_inf", infl, "_mod", toupper(substr(model, 1, 1))), fill=Ke),
       data=as.data.frame(cbind(rhat_mean10worst, sim_attrF)) #%>% filter(infl==1) 
       ) + ylab("") + geom_boxplot()

ggsave(paste0("plots/rhats10worst_", mesg, ".pdf"), height=15, width=6)



mse = sapply(resid_list, function(x) mean(abs(x^2), na.rm=TRUE))
mae = sapply(resid_list, function(x) mean(abs(x), na.rm=TRUE))
hist(mse)

mserat = sapply(mseratios_list, mean)
maerat = sapply(maeratios_list, mean)

mae_mod = lm(mae ~ model,
                 data=cbind(mae, sim_attrF))
anova(mae_mod)
summary(mae_mod)
plot(mae_mod)



ggplot(aes(x=log(mse), y=paste0("K0", K0, "_Ke", Ke, "_inf", infl, "_mod", toupper(substr(model, 1, 1))), fill=Ke),
       data=as.data.frame(cbind(mse, sim_attrF))) + ylab("") + geom_boxplot() 

ggsave(paste0("plots/mse_", mesg, ".pdf"), height=15, width=6)


ggplot(aes(x=log(mae), y=paste0("K0", K0, "_Ke", Ke, "_inf", infl, "_mod", toupper(substr(model, 1, 1))), fill=Ke),
       data=as.data.frame(cbind(mae, sim_attrF))) + ylab("") + geom_boxplot() 

ggsave(paste0("plots/mae_", mesg, ".pdf"), height=15, width=6)



ggplot(aes(x=log(mserat), y=paste0("K0", K0, "_Ke", Ke, "_inf", infl, "_mod", toupper(substr(model, 1, 1))), fill=Ke),
       data=as.data.frame(cbind(mserat, sim_attrF))) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2)

ggsave(paste0("plots/mserat_", mesg, ".pdf"), height=15, width=6)


ggplot(aes(x=log(maerat), y=paste0("K0", K0, "_Ke", Ke, "_inf", infl, "_mod", toupper(substr(model, 1, 1))), fill=Ke),
       data=as.data.frame(cbind(maerat, sim_attrF))) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2)

ggsave(paste0("plots/maerat_", mesg, ".pdf"), height=15, width=6)



ggplot(aes(x=sparse_minus_base, y=paste0("K0", K0, "_Ke", Ke, "_inf", infl), fill=Ke),
       data=pairwisediff(mse)) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2) +
  ggtitle("MSE, pairwise difference")

ggsave(paste0("plots/mse_pairwisediff_", mesg, ".pdf"), height=15, width=6)


ggplot(aes(x=sparse_minus_base, y=paste0("K0", K0, "_Ke", Ke, "_inf", infl), fill=Ke),
       data=pairwisediff(mserat)) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2) +
  ggtitle("MSE ratio, pairwise difference")

ggsave(paste0("plots/mse_ratio_pairwisediff_", mesg, ".pdf"), height=15, width=6)



ggplot(aes(x=sparse_minus_base, y=paste0("K0", K0, "_Ke", Ke, "_inf", infl), fill=Ke),
       data=pairwisediff(mae)) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2) +
  ggtitle("MAE, pairwise difference")

ggsave(paste0("plots/mae_pairwisediff_", mesg, ".pdf"), height=15, width=6)


ggplot(aes(x=sparse_minus_base, y=paste0("K0", K0, "_Ke", Ke, "_inf", infl), fill=Ke),
       data=pairwisediff(maerat)) + ylab("") + geom_boxplot() + geom_vline(xintercept=0, lty=2) +
  ggtitle("MAE ratio, pairwise difference")

ggsave(paste0("plots/mae_ratio_pairwisediff_", mesg, ".pdf"), height=15, width=6)





## Baseline profile
# library("lattice")
# i = 1
# sim_attr[i,]
# mean_prof_list[[i]]
# levelplot(log(mean_prof_list[[i]]))
# 
# 
# pow = function(x, p) x^p - 1
# pw = 0.4
# plot(mean_prof_list[[i]][1,], pow(mean_prof_list[[i]][1,], pw) + 1)


### profile plots

(elems_list = rownames(dat_stan$beta))
(src_list = unname(sapply(colnames(dat_stan$beta), function(x) strsplit(x, "_")[[1]][1])))

pow = function(x, p) x^p - 1
pw = 0.4
ticks = c(0, 0.01, 0.1, 0.25, 0.5, 1.0)

k0 = 5
pdf(file=paste0("plots/profiles_", mesg, "_k0", k0, ".pdf"), width=11, height=7)
for (inf in unique(sim_attr$infl)) {
  for (k in 1:(k0+2)) {
    for (ke in unique(sim_attr$Ke)) {
      indx_now = which(sim_attr$K0 == k0 & sim_attr$infl == inf & sim_attr$Ke == ke)
      
      if (k0 + ke >= k) {
        plot(pow(Lam_prior_list[[indx_now[1]]][,k], pw)+1, axes=F, ylab="", xlab="", ylim=c(0,1), xlim=c(1, length(elems_list)))
        axis(side=1, at=1:length(elems_list), labels = elems_list)      
        axis(side=2, at=1+pow(ticks,pw), labels=ticks)
        
        for (i in indx_now) {
          points(pow(mean_prof_list[[i]][k,], pw)+1, 
                 col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
                 pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4))
        }
        title(main=paste0(colnames(Lam_prior_list[[indx_now[1]]])[k], "\nK0=", k0, " Ke=", ke, "\ninf=", inf))
        legend("topright", pch=c(3,4,1), col=c("blue", "red", "black"), legend=c("sparse", "base", "prior"), bty="n")
      }
    }
  }
}
dev.off()



### contrib plots

ticksF = 10.0^(0:5)
nloc = 40
loc_indx = floor(seq(1, dat_stan$n, length=nloc))
# loc_indx = 1:nloc
powF = 0.2

k0 = 5
pdf(file=paste0("plots/contributions_", mesg, "_k0", k0, ".pdf"), width=11, height=7)
for (inf in unique(sim_attr$infl)) {
  for (k in 1:(k0+2)) {
    for (ke in unique(sim_attr$Ke)) {

      indx_now = which(sim_attr$K0 == k0 & sim_attr$infl == inf & sim_attr$Ke == ke)
      rhatsF_now = sapply(indx_now, function(ii) mean(rhats_list[[ii]][grep(paste0("F\\[", k, ","), names(rhats_list[[ii]]))][loc_indx]))
      
        if(k0 + ke >= k) {

        plot(NULL, axes=F, ylab="", xlab="", ylim=c(0,15), xlim=c(1,nloc))
        axis(side=1)      
        axis(side=2, at=ticksF^powF, labels=ticksF)
        
        for (i in indx_now) {
          points(med_contrib_list[[i]][loc_indx, k]^powF, 
                 col=ifelse(grepl("sparse", sim_attr[i,"model"]), "blue", "red"), 
                 pch=ifelse(grepl("sparse", sim_attr[i,"model"]), 3, 4))
        }
        
        title(main=paste0(colnames(Lam_prior_list[[indx_now[1]]])[k], "\nK0=", k0, " Ke=", ke, "\ninf=", inf))
        legend("topright", pch=c(3,4), col=c("blue", "red"), legend=c("sparse", "base"), bty="n")
        legend("topleft", ncol=3, title="Rhat", #col=rep(c("red", "blue"), times=3),
               legend=c("base", "sparse", round(tapply(rhatsF_now, sim_attr[indx_now, c("model", "Ke")], mean), 2)), bty="n")
      }
    } 
  }
}
dev.off()

