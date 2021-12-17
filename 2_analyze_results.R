# Script to analyze Bayesian model output results ##############################################################################################
################################################################################################################################################

# Load in data, do initial exploratory analysis ###############################################################################################
###############################################################################################################################################

# Load in raw Rhizoplaca lichen data file
# load(file='data/dat_rhizICP.rda')

rm(list=ls())

mod = "base"
mod = "sparse"

idf_rule = "original"

K = 5
K = 6
K = 7

inflation_factor = 1.0
inflation_factor = 2.0
inflation_factor = 3.0

seed = 12141

# load(paste0("postsim/", mod, "Model_K", K, "_idf_original_inflat_fac", inflation_factor, ".rda"))
load(paste0("postsim/", mod, "Model_idf_", idf_rule, "_K", K, "_inf", inflation_factor, "_seed", seed, ".rda"))

### notes
# base inf1 K5 idfOrig should replicate paper results (it doesn't: road dust goes all N)
# sparse inf1 K5 idfOrig should replicate paper results (it does)

###
library("rstan")

rhats = summary(fit)$summary[,"Rhat"]
length(rhats)
rhats_ord = order(rhats)
head(rhats[rhats_ord], n=50)
tail(rhats[rhats_ord], n=50)

traceplot(fit, pars='lam[22,1]')
traceplot(fit, pars='lam[4,5]')

traceplot(fit, pars='ly[66]')
traceplot(fit, pars='ly[9]')

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
plot(density(f[,1,1]))

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
dim(median_contributions) <- c(96, K)
median_contributions <- as.data.frame(median_contributions)

if (K == 5) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
} else if (K == 6) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
} else if (K == 7) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural', 'Anthropogenic')
}

median_contributions$long <- dat_now$long
median_contributions$lat <- dat_now$lat

## if base model
if (mod == "base") {
  median_contributions[,1:K] <- median_contributions[,1:K]*10000
}

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

median_contributions$long <- dat_now$long
median_contributions$lat <- dat_now$lat

if (K == 5) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust")
} else if (K == 6) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural")
} else if (K == 7) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural', 'Anthropogenic')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural", "Unspecified Anthropogenic")
}
colnames(median_contributions)[1:K] <- srcnames

k = 2

for (k in 1:K) {
  ggmap(lichen_map) +
    geom_point(data = median_contributions %>% arrange(desc(!!sym(srcnames[k]))),  # I arrange the data frame so that the larger points are plotted first (so small points don't get covered)
               aes_string(x = "long", y = "lat",
                          color = srcnames[k], size=srcnames[k])) +
    labs(color="Contributions", size="") + 
    ggtitle(srcnames_fancy[k]) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous() +
    guides(color=guide_colorbar(order=1)) +
    guides(size=guide_legend(reverse=T, order=2))
  
  pltname = paste0("plots/map_contrib_", mod, "_K", K, "_inflat", inflation_factor, "_idfOrig_", srcnames[k], ".pdf")
  ggsave(pltname, width=6.25, height=5) ## file names should reflect current selections
}

###









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
