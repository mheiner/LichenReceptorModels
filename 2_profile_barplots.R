rm(list=ls())
library("rstan")

L = 25

# idf_rule = "original"
idf_rule = "Id2" # use

K0 = 5
Ke = 2
(K = K0 + Ke) # total number, including Baseline

# inflation_factor = 0.5
inflation_factor = 1.0
# inflation_factor = 2.0
# inflation_factor = 3.0

seedB = 42702 # baseT
seedS = 42702 # sparseT

q25_prof = list()
q50_prof = list()
q75_prof = list()

mesgB = paste0("baseTModel_idf_", idf_rule, "_K0", K0, "_Ke", Ke, "_inf", inflation_factor, "_seed", seedB)
load(paste0("postsim/fit_", mesgB, ".rda"))

q25_prof[["base"]] <- summary(fit, pars='lam')$summary[,'25%']
q50_prof[["base"]] <- summary(fit, pars='lam')$summary[,'50%']
q75_prof[["base"]] <- summary(fit, pars='lam')$summary[,'75%']

dim(q25_prof[["base"]]) <- c(K, L)
dim(q50_prof[["base"]]) <- c(K, L)
dim(q75_prof[["base"]]) <- c(K, L)

mesgS = paste0("sparseTModel_idf_", idf_rule, "_K0", K0, "_Ke", Ke, "_inf", inflation_factor, "_seed", seedS)
load(paste0("postsim/fit_", mesgS, ".rda"))
elems_list = rownames(dat_stan$Lam_prior)

q25_prof[["sparse"]] <- summary(fit, pars='lam')$summary[,'25%']
q50_prof[["sparse"]] <- summary(fit, pars='lam')$summary[,'50%']
q75_prof[["sparse"]] <- summary(fit, pars='lam')$summary[,'75%']

dim(q25_prof[["sparse"]]) <- c(K, L)
dim(q50_prof[["sparse"]]) <- c(K, L)
dim(q75_prof[["sparse"]]) <- c(K, L)

colnames(q50_prof[["base"]]) <- colnames(q50_prof[["sparse"]]) <- elems_list
if(K0 == 5 & Ke == 0) {
  rownames(q50_prof[["base"]]) <- rownames(q50_prof[["sparse"]]) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
} else if (K0 == 5 & Ke == 1) {
  rownames(q50_prof[["base"]]) <- rownames(q50_prof[["sparse"]]) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
} else if (K0 == 5 & Ke == 2) {
  rownames(q50_prof[["base"]]) <- rownames(q50_prof[["sparse"]]) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural', 'Anthropogenic')
}




min_alpha = 5
source('helperScripts/all_profiles_nitrogen_modified_specificInflation_unimodal.R')
elems_list <- colnames(alpha_matrix)
(L = length(elems_list))

alpha_matrix = rbind(alpha_matrix, natural_alpha=empty_alpha, anthropogenic_alpha=empty_alpha)
beta_matrix = rbind(beta_matrix, natural_beta=natural_beta, anthropogenic_beta=anthropogenic_beta)

ls()
dim(alpha_matrix)
dim(beta_matrix)

(src_names = unname(sapply(rownames(alpha_matrix), function(x) strsplit(x, "_")[[1]][1])))
(el_names = colnames(alpha_matrix))

if(K0 == 5 & Ke == 1) {
  src_use = c(1,2,4,5,6,9)
} else if (K0 == 5 & Ke == 2) {
  src_use = c(1,2,4,5,6,9,10)
}

alpha_matrix = alpha_matrix[src_use,]
beta_matrix = beta_matrix[src_use,]

if (idf_rule == "Id2") {
  indx_zero = read.csv(file="data/hardZeroId2.csv", header=TRUE)
  nz = K - 1
  for (k in 1:K) {
    alpha_matrix[k, indx_zero[1:nz, grep( src_names[src_use][k], colnames(indx_zero))]] = 0
  }
}

source("helperScripts/simProfiles_prior.R")
dim(Lsim)

sourcenames = c("Baseline", "Playa Dust", "Brake Wear", "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural", "Unspecified Anthropogenic")

# Create a power transformation function
pow = function(x, p) x^p - 1
pw = 0.4
ticks = c(0, 0.01, 0.1, 0.25, 0.5, 1.0)
elems_list <- elem_names

k = 1
# Create and save the bar plot for each profile
for(k in 1:K) {
  prior_lower <- apply(Lsim[k,,], 1, FUN = function(x) quantile(x, 0.25))
  prior_medians <- apply(Lsim[k,,], 1, median)
  prior_upper <- apply(Lsim[k,,], 1, FUN = function(x) quantile(x, 0.75))
  
  df3 <- data.frame(Median = prior_medians, Elements = elems_list, Lower = prior_lower,
                    Upper = prior_upper, Type = "Prior")
  
  df1 <- data.frame(Median = q50_prof[["base"]][k,], Elements = elems_list,
                    Lower = q25_prof[["base"]][k,], Upper = q75_prof[["base"]][k,],
                    Type = "Base BMRM")
  
  df2 <- data.frame(Median = q50_prof[["sparse"]][k,], Elements = elems_list,
                    Lower = q25_prof[["sparse"]][k,], Upper = q75_prof[["sparse"]][k,],
                    Type = "S-BMRM")
  
  df_new <- rbind(df1, df2, df3)
  df_new$Type <- as.factor(df_new$Type)
  
  ggplot(data=df_new, aes(x=Elements, y=1+pow(Median, pw), fill=Type)) + xlab(label = "") +
    geom_bar(stat='identity', position= position_dodge())+
    geom_errorbar(aes(ymin=(1+pow(Lower, pw)), ymax=(1+pow(Upper, pw))), col='black', position = position_dodge()) +
    ggtitle(sourcenames[k]) + 
    scale_fill_discrete(type=c('red', 'gray50', 'blue'))+
    scale_y_continuous("Percent", breaks=1+pow(ticks,pw), labels=ticks*100, limits=c(0.0, 1.01))+
    theme(legend.title = element_blank())
  
  ggsave(paste0("plots/ProfileBarplots", "_K0", K0, "_Ke", Ke, "_inflat", inflation_factor, "_id", idf_rule, "_", gsub(" ", "_", sourcenames[k]), ".pdf"), 
         width=9, height=2.0)
  cat(k, "of", K, "\r")
}
