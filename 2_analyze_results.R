# Script to analyze Bayesian model output results ##############################################################################################
################################################################################################################################################

# Load in data, do initial exploratory analysis ###############################################################################################
###############################################################################################################################################

# Load in raw Rhizoplaca lichen data file
# load(file='data/dat_rhizICP.rda')

# Load in model output file being analyzed
load("postsim/baseModel_5prof_identified.rda")


# Explore estimated values
options(max.print = 7500)
print(fit, pars=c('lam[1,1]', 'lam[2,3]', 'lam[3,4]', 'lam[24,2]', 'lam[22,3]', 'lam[16,2]', 'lam[15,2]'), digits=4)
print(fit, pars='lam', digits=4)
print(fit, pars='F', digits=3)
print(fit, pars='cv', digits=3)

# Create trace plots of various estimated parameters to evaluate convergence of chains
library('rstan')
traceplot(fit, pars=c('cv[24]', 'lam[9,1]', 'lam[1,2]', 'lam[11,3]', 'lam[22,4]','F[2,50]'), inc_warmup=TRUE)
traceplot(fit, pars=c('F[1,1]', 'F[1,10]', 'F[5,15]', 'F[3,55]', 'F[2,75]', 'F[1,90]', 'F[3,6]',
                      'F[2,17]', 'F[2,18]', 'F[4,85]', 'F[5,68]', 'F[4,78]'), inc_warmup=TRUE)
traceplot(fit, pars=c('lam[1,1]', 'lam[2,3]', 'lam[3,4]', 'lam[24,2]', 'lam[13,3]', 'lam[16,5]',
                      'lam[15,5]', 'lam[25,5]'), inc_warmup=TRUE)
traceplot(fit, pars=c('cv[1]','cv[2]','cv[3]','cv[15]','cv[20]','cv[25]'), inc_warmup=TRUE)

pairs(fit, pars = c("lp__",'cv[1]', 'cv[2]', 'F[1,1]', 'lam[25,5]'), las = 1)

plot(fit, pars=c('F[2,1]', 'F[3,5]', 'F[7,10]', 'F[3,20]', 'F[5,1]', 'F[2,4]', 'F[4,7]', 'F[2,17]'))

# Look at distribution of baseline profile contributions across all MCMC samples
f <- rstan::extract(fit, pars='F')$F
plot(density(f[,1,]), xlim=c(0,25))

# Look at autocorrelation
lam155 <- rstan::extract(fit, pars='lam[15,5]')$`lam[15,5]`
acf(lam155)

###################################################################################################################################################
# Create various maps of estimated contributions ##################################################################################################
###################################################################################################################################################
library(ggmap)

load('data/dat_rhizICP.rda')

# Set up map
lichen_map <- get_stamenmap(bbox = c(left = -125, bottom = 35, right=-100, top=47), zoom=5,
                            maptype="terrain", color = "color", force = FALSE)

# Extract median contributions, create a data frame
median_contributions <- summary(fit, pars='F')$summary[ , "50%"]
dim(median_contributions) <- c(96,5)
median_contributions <- as.data.frame(median_contributions)
colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')#, "Anthropogenic")#, 'Iron_Escape')
median_contributions$long <- dat_now$long
median_contributions$lat <- dat_now$lat
#median_contributions[,1:5] <- median_contributions[,1:5]*10000

# Map out contributions with point size based on the factor of interest
ggmap(lichen_map) +
  geom_point(data = median_contributions, aes(x = long, y = lat, size = Exhaust)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() +
  ggtitle("Map of Source Contributions for Motor Vehicle Exhaust")# + 
  scale_size_area(name='Contribution',
                  max_size=8,
                  breaks=1000*c(2,4,6,8,10))
  
# Map out contributions for Baseline profile
  ggmap(lichen_map) +
    geom_point(data = median_contributions, aes(x = long, y = lat, size = Baseline)) +
    xlab("Longitude") + ylab("Latitude") + theme_bw() +
    ggtitle("Map of Source Contributions for Baseline Rhizoplaca") + 
    scale_size(name='Contribution',
               range = c(0.1, 8),
               breaks = 40000*c(1:4),
               labels = c("40000", "80000", "120000", "160000"))
  
# Clarify column names, change data frame from wide to long to allow for different visualizations
colnames(median_contributions)[1:5] <- c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                                    "Motor Vehicle Exhaust", "Unpaved Road Dust")
median_contributions$long <- dat_now$long
median_contributions$lat <- dat_now$lat
long_df <- median_contributions %>% 
  pivot_longer(c(1:5), names_to="Source", values_to="Contribution") %>%
  data.frame()

# Use facet wrap to display contributions for all sources (maps will appear as a grid)
ggmap(lichen_map) +
  geom_point(data=long_df, aes(x=long, y=lat, size=Contribution)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw() +
  ggtitle("Maps of Source Contributions") +
  scale_size_area(breaks=1000*c(1:5))+
  facet_wrap(~Source)



##########################################################################################################################
# Create new maps with both color gradient and point size based on contribution amounts ##################################

library("ggmap")
library("tidyverse")
lichen_map <- get_stamenmap(bbox = c(left = -121.5, bottom = 35.5,
                                     right=-104, top=47),
                            zoom=5, 
                            # maptype="toner-lines", 
                            maptype="terrain-background",
                            # maptype="toner-lite",
                            # maptype="terrain",
                            color = "bw", force = FALSE)

colnames(median_contributions)[1:5] <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')#, "Anthropogenic")#, 'Iron_Escape')

ggmap(lichen_map) +
  geom_point(data = median_contributions %>% arrange(desc(!!sym("Unpaved"))),  # I arrange the data frame so that the larger points are plotted first (so small points don't get covered)
             aes_string(x = "long", y = "lat",
                        color = "Unpaved", size="Unpaved")) +
  labs(color="Contributions", size="") + 
  ggtitle("Unpaved Road Dust") +
  scale_color_gradient(low = "blue", high = "red", trans='log') +
  scale_size_continuous(trans='log')+
  guides(color=guide_colorbar(order=1)) +
  guides(size=guide_legend(reverse=T, order=2))

ggsave("plots/map_sparseModel_Unpaved_Log.pdf", width=6.25, height=5) ## file names should reflect current selections


############################################################################################################################################
# Make bar plots for profile estimates #####################################################################################################
############################################################################################################################################

source("helperScripts/all_profiles_nitrogen_modified_specificInflation.R")

elems_list <- colnames(alpha_matrix)
estimated_profile_medians <- summary(fit, pars='lam')$summary[,'50%']
dim(estimated_profile_medians) <- c(5,25)
colnames(estimated_profile_medians) <- elems_list
rownames(estimated_profile_medians) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')#, "Anthropogenic")

estimated_upper <- summary(fit, pars='lam')$summary[,'75%']
dim(estimated_upper) <- c(5,25)
estimated_lower <- summary(fit, pars='lam')$summary[,'25%']
dim(estimated_lower) <- c(5,25)

## save these for each model run, e.g., 
# save("postsim/base_profile_medians_intervals.rda", 
#      estimated_profile_medians, estimated_profile_lower, estimated_profile_upper)


#######################################################################################################################
# Create prior vs posterior profile bar plots (with uncertainty) with re-scaled y-axis ################################

# First, load in relevant data (if not done already)
load('postsim/base_profile_medians_intervals.rda') # saved medians, upper, and lower (from commented code above)
estimated_medians_base <- estimated_profile_medians
estimated_upper_base <- estimated_upper; estimated_lower_base <- estimated_lower

load('postsim/sparse_profile_medians_intervals.rda') # saved medians, upper, and lower (from commented code above)
estimated_medians_sparse <- estimated_profile_medians
estimated_upper_sparse <- estimated_upper; estimated_lower_sparse <- estimated_lower

source('helperScripts/sim_profile_prior.R')

sourcenames = c("Baseline", "Playa Dust", "Brake Wear", "Motor Vehicle Exhaust", "Unpaved Road Dust")
simnames = c("baseline_sims", "playa_sims", "brake_sims", "exhaust_sims", "unpaved_dust_sims")

# Create a power transformation function
pow = function(x, p) x^p - 1
pw = 0.4
ticks = c(0, 0.01, 0.1, 0.25, 0.5, 1.0)
elems_list <- elem_names

profile = 1

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
  
  df_new <- rbind(df1, df2, df3)
  df_new$Type <- as.factor(df_new$Type)
  
  ggplot(data=df_new, aes(x=Elements, y=1+pow(Median, pw), fill=Type)) + xlab(label = "") +
    geom_bar(stat='identity', position= position_dodge())+
    geom_errorbar(aes(ymin=(1+pow(Lower, pw)), ymax=(1+pow(Upper, pw))), col='black', position = position_dodge()) +
    ggtitle(sourcenames[profile]) + 
    scale_fill_discrete(type=c('red', 'gray50', 'blue'))+
    scale_y_continuous("Percent", breaks=1+pow(ticks,pw), labels=ticks*100, limits=c(0.0, 0.96))+
    theme(legend.title = element_blank())
  
  ggsave(paste0("plots/ProfileBarplots_", gsub(" ", "_", sourcenames[profile]), ".pdf"), width=9, height=2.0)
  cat(profile, "of", length(sourcenames), "\r")
}





#######################################################################################################################
# Analyze other fitted values (residuals, total contributions, etc.) ##################################################
#######################################################################################################################

same_elems <- colnames(dat_now) %in% colnames(estimated_profile_medians)

estimated_y <- exp(summary(fit, pars='lLF')$summary[, "50%"]) * 10000 # scale estimates up to ppm if using base model
dim(estimated_y) <- c(96,25)
raw_data = dat_now[,elems_list]

colnames(estimated_y) <- colnames(raw_data)

resid_mat <- raw_data - estimated_y

resid_means <- apply(resid_mat, 2, FUN=function(x) mean(x, na.rm=T))
resid_sds <- apply(resid_mat, 2, FUN=function(x) sd(x, na.rm=T))

elem = 1
hist(resid_mat[,elem], xlab="Residuals",
     main=paste0("Histogram of Residuals for ", elems_list[elem]))

hist( ((resid_mat[,elem] - resid_means[elem]) / resid_sds[elem] ), xlab="Standardized Residuals",
      main=paste0("Histogram of Standardized Residuals for ", elems_list[elem]))


# Plot fitted vs observed values for each element
elem=1
df <- data.frame(log_est = log(estimated_y[,elem]), log_raw = log(raw_data[,elem]))
ggplot(data=df, aes(x=log_est, y=log_raw)) + geom_point() +
  xlab(paste0("log of estimated ", elems_list[elem])) + ylab(paste0('log of true ', elems_list[elem])) + geom_smooth(method='lm',se=F)


lresid <- log(raw_data) - log(estimated_y)

resid_ind <- which(!is.na(lresid[,elem]))
raw_ind <- which(!is.na(raw_data[,elem]))

if(all(resid_ind==raw_ind)){
  1 - ( sum(lresid[resid_ind,elem]^2) / sum( (log(raw_data[raw_ind,elem]) - mean(log(raw_data[raw_ind,elem])))^2 ) )
}else{
  print("Nope")
}

# Compute a measure of r^2 for each element
model_r2 <- numeric(length=25)
for(i in 1:25){
  elem=i
  
  lresid <- log(raw_data) - log(estimated_y)
  
  resid_ind <- which(!is.na(lresid[,elem]))
  raw_ind <- which(!is.na(raw_data[,elem]))
  
  if(all(resid_ind==raw_ind)){
    model_r2[i] <- 1 - ( sum(lresid[resid_ind,elem]^2) / sum( (log(raw_data[raw_ind,elem]) - mean(log(raw_data[raw_ind,elem])))^2 ) )
  }else{
    print("Nope")
  }
  
}
names(model_r2) <- elems_list
round(model_r2[sort(elems_list)],2)

# Total mass (estimate / true value)
mass_explained <- ( (1/exp(lresid))*100 ) - 100
ggplot() + geom_histogram(aes(x=mass_explained[abs(mass_explained)<100])) + ylab("Count") +
  xlab("(Estimated Mass / True Mass) for each Element") +
  ggtitle("Histogram of the Ratio of Estimated Mass to True Mass for each Element")

#################################################################################################################
# Analyze estimated total contributions #########################################################################

median_contributions <- median_contributions[,1:5]
ggplot() + geom_histogram(aes(x=rowSums(median_contributions)), bins=8) + xlab("Total Estimated Mass per Sample") +
  ggtitle("Histogram of Total Estimated Mass for each Sample") + ylab("Count")


## save these for each model run, e.g., 
# save("postsim/base_contribution_medians.rda", 
#      estimated_contributions)

load('postsim/base_contribution_medians.rda') # saved medians (from commented code above)
median_contrib_base <- median_contributions

load('postsim/sparse_profile_medians_intervals.rda') # saved medians (from commented code above)
median_contrib_sparse <- median_contributions


# Compare Base vs Sparse sample-wise total contributions (histogram)
ggplot() +
  geom_histogram(mapping = aes(x = rowSums(median_contrib_base)/rowSums(median_contrib_sparse)))+
  ggtitle("Ratio of Sample-wise Base vs Sparse Total Contributions")+
  xlab("(Base / Sparse) Sample-wise Total Contributions")+
  ylab("Count")

# Also look at the same thing as above, but excluding baseline
ggplot() +
  geom_histogram(mapping = aes(x = rowSums(median_contrib_base[,-1])/rowSums(median_contrib_sparse[,-1])))+
  ggtitle("Ratio of Sample-wise Base vs Sparse Total Contributions (excluding Baseline)")+
  xlab("(Base / Sparse) Sample-wise Total Contributions (excluding Baseline)")+
  ylab("Count")+
  xlim(0,30)

# Just baseline
sum(median_contrib_base[,1]) / sum(median_contrib_sparse[,1])
