rm(list=ls())

idf_rule = "original"
K = 6
inflation_factor = 2.0

source("helperScripts/all_profiles_nitrogen_modified_specificInflation.R")
elems_list <- colnames(alpha_matrix)
(L = length(elems_list))


load(paste0("postsim/baseModel_K", K, "_idf_original_inflat_fac", inflation_factor, ".rda"))
load(paste0("postsim/baseModel_idf_", idf_rule, "_K", K, "_inf", inflation_factor, "_seed", seed, ".rda"))


estimated_medians_base <- summary(fit, pars='lam')$summary[,'50%']
dim(estimated_medians_base) <- c(K,L)
estimated_upper_base <- summary(fit, pars='lam')$summary[,'75%']
dim(estimated_upper_base) <- c(K,L)
estimated_lower_base <- summary(fit, pars='lam')$summary[,'25%']
dim(estimated_lower_base) <- c(K,L)


load(paste0("postsim/sparseModel_K", K, "_idf_original_inflat_fac", inflation_factor, ".rda"))
load(paste0("postsim/sparseModel_idf_", idf_rule, "_K", K, "_inf", inflation_factor, "_seed", seed, ".rda"))

estimated_medians_sparse <- summary(fit, pars='lam')$summary[,'50%']
dim(estimated_medians_sparse) <- c(K,L)
estimated_upper_sparse <- summary(fit, pars='lam')$summary[,'75%']
dim(estimated_upper_sparse) <- c(K,L)
estimated_lower_sparse <- summary(fit, pars='lam')$summary[,'25%']
dim(estimated_lower_sparse) <- c(K,L)


colnames(estimated_medians_base) <- colnames(estimated_medians_sparse) <- elems_list
if(K == 5) {
  rownames(estimated_medians_base) <- rownames(estimated_medians_sparse) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
} else if (K == 6) {
  rownames(estimated_medians_base) <- rownames(estimated_medians_sparse) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
}


source('helperScripts/sim_profile_prior.R')

if(K == 5) {
  sourcenames = c("Baseline", "Playa Dust", "Brake Wear", "Motor Vehicle Exhaust", "Unpaved Road Dust")
  simnames = c("baseline_sims", "playa_sims", "brake_sims", "exhaust_sims", "unpaved_dust_sims")
} else if (K == 6) {
  sourcenames = c("Baseline", "Playa Dust", "Brake Wear", "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural")
  simnames = c("baseline_sims", "playa_sims", "brake_sims", "exhaust_sims", "unpaved_dust_sims", "natural_sims")
}


# Create a power transformation function
pow = function(x, p) x^p - 1
pw = 0.4
ticks = c(0, 0.01, 0.1, 0.25, 0.5, 1.0)
elems_list <- elem_names

profile = 5

# Create and save the bar plot for each profile
for(profile in 1:length(sourcenames)) {
  sims_now <- get(simnames[profile])
  
  prior_lower <- apply(sims_now, 1, FUN = function(x) quantile(x, 0.25))
  prior_medians <- apply(sims_now, 1, median)
  prior_upper <- apply(sims_now, 1, FUN = function(x) quantile(x, 0.75))
  
  df3 <- data.frame(Median = prior_medians, Elements = elems_list, Lower = prior_lower,
                    Upper = prior_upper, Type = "Prior")
  
  df1 <- data.frame(Median = estimated_medians_base[profile,], Elements = elems_list,
                    Lower = estimated_lower_base[profile,], Upper = estimated_upper_base[profile,],
                    Type = "Base BMRM")
  
  df2 <- data.frame(Median = estimated_medians_sparse[profile,], Elements = elems_list,
                    Lower = estimated_lower_sparse[profile,], Upper = estimated_upper_sparse[profile,],
                    Type = "S-BMRM")
  
  ## pick one
  df_new <- rbind(df1, df2, df3)
  df_new$Type <- as.factor(df_new$Type)
  
  # df_new <- rbind(df1, df3)
  # df_new$Type <- as.factor(df_new$Type)
  
  # df_new <- rbind(df2, df3)
  # df_new$Type <- as.factor(df_new$Type)
  
  ggplot(data=df_new, aes(x=Elements, y=1+pow(Median, pw), fill=Type)) + xlab(label = "") +
    geom_bar(stat='identity', position= position_dodge())+
    geom_errorbar(aes(ymin=(1+pow(Lower, pw)), ymax=(1+pow(Upper, pw))), col='black', position = position_dodge()) +
    ggtitle(sourcenames[profile]) + 
    scale_fill_discrete(type=c('red', 'gray50', 'blue'))+
    scale_y_continuous("Percent", breaks=1+pow(ticks,pw), labels=ticks*100, limits=c(0.0, 1.01))+
    theme(legend.title = element_blank())
  
  ggsave(paste0("plots/ProfileBarplots", "_K", K, "_inflat", inflation_factor, "_idfOrig", "_", gsub(" ", "_", sourcenames[profile]), ".pdf"), 
         width=9, height=2.0)
  cat(profile, "of", length(sourcenames), "\r")
}

