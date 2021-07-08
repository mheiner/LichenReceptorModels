rm(list=ls())

library(tidyverse)
require(maps)
require(ggmap)
require(viridis)

source("helperScripts/detection_limits.R")
dat_now0 = read.csv("data/rhizoplaca96.csv", header=TRUE)
head(dat_now0)
tail(dat_now0)

el_indx = c(7:35, 38) # indexes of elemental concentrations
perc_elements = c("Ca", "K", "Mg", "N", "P", "S")

## BDL, missing values
dat_now0_missing = dat_now0[,el_indx] %>% is.na() | dat_now0[,el_indx] == ""
dat_now0_nd = apply(dat_now0[,el_indx], c(1,2), function(x) grepl("nd", x))
dat_now0_bdl = apply(dat_now0[,el_indx], c(1,2), function(x) grepl("bdl", x))

colMeans(dat_now0_missing)
colMeans(dat_now0_nd)
colMeans(dat_now0_bdl)

(propmiss = colMeans(dat_now0_missing + dat_now0_nd + dat_now0_bdl))
dat_now0[,el_indx[propmiss > 0.1]]
dat_now0[,el_indx[propmiss > 0.8]]

el_use_indx = el_indx[which(propmiss <= 0.8)]

colnames(dat_now0)[el_indx]
colnames(dat_now0)[el_use_indx]
setdiff(colnames(dat_now0)[el_indx], colnames(dat_now0)[el_use_indx]) # excluded elements

length(el_use_indx)




## PERC to PPM
dat_now = dat_now0
dat_now[,el_indx] = apply(dat_now[,el_indx], c(1,2), function(x) gsub("[a-z]", NA, x)) ## remove letters
dat_now[,el_indx] = apply(dat_now[,el_indx], 2, as.numeric)
head(dat_now)
str(dat_now)

dat_now[,grep("PERC", colnames(dat_now))] = 1.0e4 * dat_now[,grep("PERC", colnames(dat_now))]
colnames(dat_now) = gsub("PERC", "", colnames(dat_now))
colnames(dat_now0_missing) = gsub("PERC", "", colnames(dat_now0_missing))
colnames(dat_now0_nd) = gsub("PERC", "", colnames(dat_now0_nd))
colnames(dat_now0_bdl) = gsub("PERC", "", colnames(dat_now0_bdl))

head(dat_now)
(el_use = colnames(dat_now)[el_use_indx])

## Detection limits
DLtab = cbind(min=apply(dat_now[,el_indx], 2, min, na.rm=TRUE)[names(detection_limits_ppmDilute)], DL=detection_limits_ppmDilute)
cbind(DLtab, diff=DLtab[,1] - DLtab[,2], ratio=round(DLtab[,1] / DLtab[,2], 4))

(DL = sapply( colnames(dat_now[el_use_indx]), function(elem) { 
  ifelse(elem %in% names(detection_limits_ppmDilute), detection_limits_ppmDilute[elem], min(dat_now[,elem], na.rm=TRUE)) 
  } ))



## indexing 
dim(dat_now[,el_use_indx])
dim(dat_now0_missing[,el_use])

colSums(dat_now0_missing[,el_use])
colSums(dat_now0_nd[,el_use])
colSums(dat_now0_bdl[,el_use])


Y = as.matrix(dat_now[,el_use_indx])
(n = nrow(Y))
(L = ncol(Y))
n*L


### Perform indexing on the L by n transpose of Y

indx_obs = which(!is.na(t(Y)))
indx_mis = which(t(dat_now0_missing[,el_use]))
indx_bdl = which(t(dat_now0_bdl[,el_use]))
dl_censor_vals = (t(dat_now0_bdl[,el_use]) * t(matrix(rep(DL[el_use], each=n), nrow=n)) )[indx_bdl]

(n_bdl = length(indx_bdl))
(n_mis = length(indx_mis))
(n_obs = length(indx_obs))

stopifnot( (n_bdl + n_mis + n_obs) == n*L )
y_obs = t(Y)[indx_obs]



dat_stan = list(n=n, L=L, K=NULL,
                n_obs=n_obs, indx_obs=indx_obs, y_obs=y_obs,
                n_mis=n_mis, indx_mis=indx_mis,
                n_bdl=n_bdl, indx_bdl=indx_bdl, dl_censor_vals=dl_censor_vals,
                alpha=NULL, beta=NULL,
                a=NULL, b=NULL,
                flim=175000,
                gamma_0=30000,
                note="Indexing is done on t(Y), which is L by n")

save(file="data/dat_rhizICP.rda", dat_stan, dat_now)
