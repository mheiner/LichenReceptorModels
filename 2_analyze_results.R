## This script is intended for use to explore results interactively, NOT sequentially or in batch mode

rm(list=ls())
source("0_prep_data.R")
(N = nrow(Y))
dat_now$ycollect[which(dat_now$ycollect == 1992)] = NA

mod = "base"
mod = "baseT"
mod = "sparse"
mod = "sparseT"

# idf_rule = "original"
idf_rule = "Id2" # use

K0 = 5
Ke = 2
(K = K0 + Ke) # total number of sources, including Baseline

# inflation_factor = 0.5
inflation_factor = 1.0
# inflation_factor = 2.0
# inflation_factor = 3.0

seed = 42702 # baseT
seed = 42702 # sparseT

mesg = paste0(mod, "Model_idf_", idf_rule, "_K0", K0, "_Ke", Ke, "_inf", inflation_factor, "_seed", seed)


load(paste0("postsim/fit_", mesg, ".rda"))
library("rstan")

(src_names = sapply(strsplit(colnames(dat_stan$beta), "_"), function(xx) xx[1]))
(elem_names = rownames(dat_stan$beta))

median_contributions <- summary(fit, pars='F')$summary[ , "50%"]
mean_profiles <- summary(fit, pars='lam')$summary[,'mean']

head(median_contributions, n=30)
head(mean_profiles, n=30)

dim(mean_profiles) <- c(K, L)
dim(median_contributions) <- c(N, K)


lFL_hat = matrix(summary(fit, pars="lLF")$summary[,"50%"], nrow=N) # estimate of log median elemental concentrations
FL_hat = exp(lFL_hat)

dim(lFL_hat)
dim(Y)

lE = log(Y) - lFL_hat # log residuals
lE0 = log(Y) - log( tcrossprod( rep(1, N), colMeans(Y, na.rm=TRUE) ) ) # log residuals from a model using only column means of Y

(mseratios = colMeans(lE^2, na.rm=TRUE) / colMeans(lE0^2, na.rm=TRUE))
plot(mseratios); abline(h=1, lty=2)

(maeratios = colMeans(abs(lE), na.rm=TRUE) / colMeans(abs(lE0), na.rm=TRUE))
plot(maeratios); abline(h=1, lty=2)


## spatial and temporal variation in residuals, contributions
datCE = dat_now %>% bind_cols(rename_with(data.frame(lE), function(x) paste0("lE_", x))) %>% 
  bind_cols(rename_with(data.frame(median_contributions), function(x) paste(src_names)))
head(datCE)


k = 2
for ( k in 1:K ) {
  plot(datCE$ycollect, datCE[,src_names[k]], main=src_names[k]) # ycollect is year collected, source refers to estimated contributions
}

for ( k in 1:K ) {
  boxplot( log(datCE[,src_names[k]]) ~ datCE$ycollect, main=src_names[k])
}

for ( ell in 1:L ) { # ycollect is confounded w location
  plot(datCE$ycollect, datCE[, paste0("lE_", elem_names[ell])], main=paste0("log resids: ", elem_names[ell]))
  abline(h=0, lty=2)
}

ell = 22
for ( ell in 1:L ) {
  plot(datCE$long, datCE[, paste0("lE_", elem_names[ell])], main=paste0("log resids: ", elem_names[ell]))
  abline(h=0, lty=2)
}




### check MCMC
library("rstan")

rhats = summary(fit)$summary[,"Rhat"]
length(rhats)
rhats_ord = order(rhats)
head(rhats[rhats_ord], n=50)
tail(rhats[rhats_ord], n=50)

traceplot(fit, pars='lam[1,2]')
traceplot(fit, pars='lam[7,4]')
traceplot(fit, pars='lLF[24,31]')
traceplot(fit, pars="F[4,77]")
traceplot(fit, pars="cv[24]")

# traceplot(fit, pars='ly[66]')
# traceplot(fit, pars='ly[9]')

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
traceplot(fit, pars="nu")

pairs(fit, pars = c("lp__",'cv[1]', 'cv[2]', 'F[1,1]', 'lam[25,5]'), las = 1)

plot(fit, pars=c('F[2,1]', 'F[3,5]', 'F[5,10]', 'F[3,20]', 'F[5,1]', 'F[2,4]', 'F[4,7]', 'F[2,17]'))
plot(fit, pars=c('F[1,1]', 'F[1,5]', 'F[1,10]', 'F[1,20]', 'F[1,32]', 'F[1,40]', 'F[1,59]', 'F[1,60]'))
plot(fit, pars=c('F[2,11]', 'F[2,51]', 'F[2,23]', 'F[2,21]', 'F[2,33]', 'F[2,44]', 'F[2,55]', 'F[2,66]'))
plot(fit, pars='cv')

# Look at distribution of baseline profile contributions across all MCMC samples
f <- rstan::extract(fit, pars='F')$F
plot(density(f[,1,1]))

# Look at autocorrelation
lam155 <- rstan::extract(fit, pars='lam[15,5]')$`lam[15,5]`
acf(lam155)




###################################################################################################################################################
# Create various maps of estimated contributions ##################################################################################################
###################################################################################################################################################

library("ggmap")
load('data/dat_rhizICP.rda')

# Set up map
lichen_map <- get_stamenmap(bbox = c(left = -125, bottom = 35, right=-100, top=47), zoom=5,
                            maptype="terrain", color = "color", force = FALSE)

# Extract median contributions, create a data frame
median_contributions <- summary(fit, pars='F')$summary[ , "50%"]
dim(median_contributions) <- c(N, K)
median_contributions <- as.data.frame(median_contributions)

if (K0 == 5 & Ke == 0) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
} else if (K0 == 5 & Ke == 1) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
} else if (K0 == 5 & Ke == 2) {
  colnames(median_contributions) <- c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural', 'Anthropogenic')
}

median_contributions$long <- dat_now$long
median_contributions$lat <- dat_now$lat

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
median_contributions$year <- dat_now$ycollect

if (K0 == 5 & Ke == 0) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust")
} else if (K0 == 5 & Ke == 1) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural")
} else if (K0 == 5 & Ke == 2) {
  srcnames = c('Baseline', 'Playa', 'Brake', 'Exhaust', 'Unpaved', 'Natural', 'Anthropogenic')
  srcnames_fancy = c("Baseline Rhizoplaca", "Playa", "Brake Wear",
                     "Motor Vehicle Exhaust", "Unpaved Road Dust", "Unspecified Natural", "Unspecified Anthropogenic")
}
colnames(median_contributions)[1:K] <- srcnames


k = 5

for (k in 1:K) {
  ggmap(lichen_map) +
    geom_point(data = median_contributions %>% arrange(desc(!!sym(srcnames[k]))),  # arrange the data frame so that the larger points are plotted first (so small points don't get covered)
               aes_string(x = "long", y = "lat",
                          color = srcnames[k], size=srcnames[k])) +
    labs(color="Contributions", size="") + 
    ggtitle(srcnames_fancy[k]) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous() +
    guides(color=guide_colorbar(order=1)) +
    guides(size=guide_legend(reverse=T, order=2))
  
  pltname = paste0("plots/map_contrib_", mesg, "_src_", srcnames[k], ".pdf")
  ggsave(pltname, width=6.25, height=5) ## file names should reflect current selections
}

### concentrations and residuals by year

ggmap(lichen_map) +
  geom_point(data = datCE %>% arrange(desc(ycollect)) %>% mutate(ycollect=as.factor(ycollect)),  # arrange the data frame so that the larger points are plotted first (so small points don't get covered)
             aes(x=long, y=lat, color = ycollect), size=3) + labs(color="year\ncollected")
ggsave("plots/yearMap.pdf", width=6, height=5)


k = 1
for (k in 1:K) {
  pp = ggmap(lichen_map) +
    geom_point(data = datCE %>% arrange(desc(!!sym(src_names[k]))),  # arrange the data frame so that the larger points are plotted first (so small points don't get covered)
               aes_string(x = "long", y = "lat",
                          color = "ycollect", size=src_names[k])) +
    labs(color="Year", size="") + 
    ggtitle(srcnames_fancy[k]) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous() +
    guides(color=guide_colorbar(order=1)) +
    guides(size=guide_legend(reverse=T, order=2))
  print(pp)
}

ell = 1
for (ell in 1:L) {
  pp = ggmap(lichen_map) +
    geom_point(data = datCE %>% arrange(desc(!!sym(paste0("lE_", elem_names[ell])))),  # arrange the data frame so that the larger points are plotted first (so small points don't get covered)
               aes_string(x = "long", y = "lat",
                          color = "ycollect", size=paste0("lE_", elem_names[ell]) )) +
    labs(color="Year", size="") + 
    ggtitle(paste0("log resid ", elem_names[ell])) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_binned(n.breaks=10) +
    guides(color=guide_colorbar(order=1)) +
    guides(size=guide_legend(reverse=T, order=2))
  print(pp)
}


### spatial variation in concentrations, residuals
library("gstat")
library("sp")

datCEsp = datCE

k = 7
for (k in 1:K) {
  datCEsp = datCE %>% select(long, lat, sym(src_names[k]))
  vg = variogram( as.formula(paste0(src_names[k], " ~ 1")), data=datCEsp, 
                  locations = ~ long + lat)
  print(plot(vg, main=src_names[k]))
}

# ell = 1
# for (ell in 1:L) {
#   datCEsp = datCE %>% drop_na(paste0("lE_", elem_names[ell]))
#   vg = variogram( as.formula(paste0("lE_", elem_names[ell], " ~ 1")), 
#                   data=datCEsp, locations = ~ long + lat)
#   print(plot(vg, main=elem_names[ell]))
# }

library("gridExtra")
ell = 1

for (ell in 1:L) { # log residual variograms and maps for each element
  
  # remove outliers
  zz0 = datCE[,paste0("lE_", elem_names[ell])]
  zzc = scale(zz0)[,1]
  datCEsp = datCE[which(abs(zzc) < 2.5),]
  
  pp = ggmap(lichen_map) +
    geom_point(data = datCEsp %>% arrange(desc(!!sym(paste0("lE_", elem_names[ell])))),  # arrange the data frame so that the larger points are plotted first (so small points don't get covered)
               aes_string(x = "long", y = "lat",
                          color = paste0("lE_", elem_names[ell])), size=1.5) +
    labs(color="log resid") + 
    ggtitle(paste0("log residuals ", elem_names[ell])) +
    scale_color_distiller(palette="Spectral") +
    guides(color=guide_colorbar(order=1)) +
    guides(size=guide_legend(reverse=T, order=2))
  
  datCEsp = datCEsp %>% drop_na(paste0("lE_", elem_names[ell]))
  vg = variogram( as.formula(paste0("lE_", elem_names[ell], " ~ 1")), 
                  data=datCEsp, locations = ~ long + lat)
  # plot(vg, main=elem_names[ell])
  
  ppvg = ggplot(vg, aes(x=dist, y=gamma)) + geom_point() + 
    ggtitle(elem_names[ell]) + ylim(c(0, max(vg$gamma))) +
    ylab("semivariance") + xlab("distance (degrees)")
  
  gg = grid.arrange(pp, ppvg, ncol=2, layout_matrix=matrix(c(1,1,1, 1,1,1, NA,2,NA, NA,2,NA), nrow=3)  )  

  ggsave(gg, file=paste0("plots/residMap_", elem_names[ell], "_", mesg, ".pdf"), width=8, height=5)
  
  rm(zz0, zzc, pp, ppvg, vg, datCEsp)
}
## signals in Ti, Na, Mo, Mn, Fe, Cu, Cr, Cd, Ba, B, As, Al, S, K




#######################################################################################################################
# Analyze other fitted values (residuals, total contributions, etc.) ##################################################
#######################################################################################################################

# same_elems <- colnames(dat_now) %in% colnames(mean_profiles)
elems_list = colnames(Y)

estimated_y = FL_hat
colnames(estimated_y) <- elems_list

resid_mat <- Y - estimated_y

resid_means <- apply(resid_mat, 2, FUN=function(x) mean(x, na.rm=T))
resid_sds <- apply(resid_mat, 2, FUN=function(x) sd(x, na.rm=T))

elem = 2
hist(resid_mat[,elem], xlab="Residuals",
     main=paste0("Histogram of Residuals for ", elems_list[elem]))

for(ell in 1:L) hist(lE[,ell], main=elems_list[ell])
for(ell in 1:L) {
  qqnorm(lE[,ell], main=elems_list[ell]); qqline(lE[,ell])
}

hist( ((resid_mat[,elem] - resid_means[elem]) / resid_sds[elem] ), xlab="Standardized Residuals",
      main=paste0("Histogram of Standardized Residuals for ", elems_list[elem]))


# Plot fitted vs observed values for each element
elem=7
for (elem in 1:L) {
  df <- data.frame(log_est = log(estimated_y[,elem]), log_raw = log(Y[,elem]))
  pp = ggplot(data=df, aes(x=log_est, y=log_raw)) + geom_point() +
    xlab(paste0("log of estimated ", elems_list[elem])) + ylab(paste0('log of true ', elems_list[elem])) + geom_smooth(method='lm', se=F)
  print(pp)
}


lresid <- log(Y) - log(estimated_y)

resid_ind <- which(!is.na(lresid[,elem]))
raw_ind <- which(!is.na(Y[,elem]))

if(all(resid_ind==raw_ind)){
  1 - ( sum(lresid[resid_ind,elem]^2) / sum( (log(Y[raw_ind,elem]) - mean(log(Y[raw_ind,elem])))^2 ) )
}else{
  print("Nope")
}


# Compute a measure of r^2 for each element
model_r2 <- numeric(length=L)
for(i in 1:L){
  elem=i

  resid_ind <- which(!is.na(lresid[,elem]))
  raw_ind <- which(!is.na(Y[,elem]))
  
  if(all(resid_ind==raw_ind)){
    model_r2[i] <- 1 - ( sum(lresid[resid_ind,elem]^2) / sum( (log(Y[raw_ind,elem]) - mean(log(Y[raw_ind,elem])))^2 ) )
  }else{
    print("Nope")
  }
  
}
round(1-mseratios,2)

names(model_r2) <- elems_list
round(model_r2[sort(elems_list)],2)
round(1-mseratios[sort(elems_list)],2)



## compare estimated mass between base and sparse models

library("rstan")

median_contrib = list()
mean_prof = list()
mean_cv = list()

seedB = 42702 # baseT
seedS = 42702 # sparseT


mesgB = paste0("baseTModel_idf_", idf_rule, "_K0", K0, "_Ke", Ke, "_inf", inflation_factor, "_seed", seedB)
load(paste0("postsim/fit_", mesgB, ".rda"))

median_contrib[["base"]] <- summary(fit, pars='F')$summary[ , "50%"] * 10e3
mean_prof[["base"]] <- summary(fit, pars='lam')$summary[,'mean']

dim(mean_prof[["base"]]) <- c(K, L)
dim(median_contrib[["base"]]) <- c(N, K)

mean_cv[["base"]] <- summary(fit, pars='cv')$summary[,'mean']


mesgS = paste0("sparseTModel_idf_", idf_rule, "_K0", K0, "_Ke", Ke, "_inf", inflation_factor, "_seed", seedS)
load(paste0("postsim/fit_", mesgS, ".rda"))

median_contrib[["sparse"]] <- summary(fit, pars='F')$summary[ , "50%"]
mean_prof[["sparse"]] <- summary(fit, pars='lam')$summary[,'mean']

dim(mean_prof[["sparse"]]) <- c(K, L)
dim(median_contrib[["sparse"]]) <- c(N, K)

mean_cv[["sparse"]] <- summary(fit, pars='cv')$summary[,'mean']


# Total mass (estimate / true value)
mass_explained <- ( (1/exp(lresid))*100 ) - 100 # 1/exp(lresid) is equal to Yhat/Y
ggplot() + geom_histogram(aes(x=mass_explained[abs(mass_explained)<100])) + ylab("Count") +
  xlab("(Estimated Mass / True Mass) for each Element") +
  ggtitle("Histogram of the Ratio of Estimated Mass to True Mass for each Element")

#################################################################################################################
# Analyze estimated total contributions #########################################################################

# median_contributions <- median_contributions[,1:5]
ggplot() + geom_histogram(aes(x=rowSums(median_contributions)), bins=8) + xlab("Total Estimated Mass per Sample") +
  ggtitle("Histogram of Total Estimated Mass for each Sample") + ylab("Count")


# Compare Base vs Sparse sample-wise total contributions (histogram)
ggplot() +
  geom_histogram(mapping = aes(x = rowSums(median_contrib[["base"]])/rowSums(median_contrib[["sparse"]]))) +
  geom_vline(xintercept=1, lty=2) +
  ggtitle("Ratio of Sample-wise Base vs Sparse Total Contributions") +
  xlab("(Base / Sparse) Sample-wise Total Contributions") +
  ylab("Count") 

# Also look at the same thing as above, but excluding baseline
ggplot() +
  geom_histogram(mapping = aes(x = rowSums(median_contrib[["base"]][,-1])/rowSums(median_contrib[["sparse"]][,-1])))+
  geom_vline(xintercept=1, lty=2) +
  ggtitle("Ratio of Sample-wise Base vs Sparse Total Contributions (excluding Baseline)")+
  xlab("(Base / Sparse) Sample-wise Total Contributions (excluding Baseline)")+
  ylab("Count")+
  xlim(0,30)

mean( rowSums(median_contrib[["base"]][,-1])/rowSums(median_contrib[["sparse"]][,-1]) > 10 )

# Just baseline
sum(median_contrib[["base"]][,1]) / sum(median_contrib[["sparse"]][,1])

plot(mean_cv[["base"]], mean_cv[["sparse"]], pch=20)
abline(0, 1, lty=2)
text(mean_cv[["base"]], mean_cv[["sparse"]]+0.02, labels=elems_list)
mean(mean_cv[["base"]] > mean_cv[["sparse"]])
