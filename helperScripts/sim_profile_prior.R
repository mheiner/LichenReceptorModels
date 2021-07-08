
# Create a function to simulate means & sd's for each profile (for a single iteration)
single_draw <- function(alpha_vector, beta_vector){
  K <- length(alpha_vector)
  stopifnot(length(beta_vector) == K)
  xx <- rgamma(K, shape = alpha_vector, rate = beta_vector)

  lambda_out <- xx / sum(xx)
  
  return(lambda_out)
}

set.seed(5)

#####################################################################################################################################
# Now, go through each profile, getting 100,000 draws for each profile, so that normalized medians, means, and sds can be found #####
#####################################################################################################################################


## Must first source '0_all_profiles_nitrogen_modified_specificInflation.R'
# source("all_profiles_nitrogen_modified_specificInflation.R")


# Agriculture Soil
agriculture_soil_sims <- replicate(1e5, single_draw(agriculture_alpha, agriculture_beta))

agriculture_sim_medians <- apply(agriculture_soil_sims, 1, FUN = median)
agriculture_sim_means <- apply(agriculture_soil_sims, 1, FUN = mean)
agriculture_sim_sds <- apply(agriculture_soil_sims, 1, FUN = sd)

# Coal Combustion
coal_combustion_sims <- replicate(1e5, single_draw(c_c_alpha, c_c_beta))

coal_combustion_sim_medians <- apply(coal_combustion_sims, 1, FUN = median)
coal_combustion_sim_means <- apply(coal_combustion_sims, 1, FUN = mean)
coal_combustion_sim_sds <- apply(coal_combustion_sims, 1, FUN = sd)

# Unpaved Road Dust
unpaved_dust_sims <- replicate(1e5, single_draw(unpaved_alpha, unpaved_beta))

unpaved_dust_sim_medians <- apply(unpaved_dust_sims, 1, FUN = median)
unpaved_dust_sim_means <- apply(unpaved_dust_sims, 1, FUN = mean)
unpaved_dust_sim_sds <- apply(unpaved_dust_sims, 1, FUN = sd)

# Coal-Fired Power Plant
coal_plant_sims <- replicate(1e5, single_draw(coal_plant_alpha, coal_plant_beta))

coal_plant_sim_medians <- apply(coal_plant_sims, 1, FUN = median)
coal_plant_sim_means <- apply(coal_plant_sims, 1, FUN = mean)
coal_plant_sim_sds <- apply(coal_plant_sims, 1, FUN = sd)

# Oil Refinery
oil_refinery_sims <- replicate(1e5, single_draw(oil_alpha, oil_beta))

oil_refinery_sim_medians <- apply(oil_refinery_sims, 1, FUN = median)
oil_refinery_sim_means <- apply(oil_refinery_sims, 1, FUN = mean)
oil_refinery_sim_sds <- apply(oil_refinery_sims, 1, FUN = sd)

# Avg. Fine Playa

playa_sims <- replicate(1e5, single_draw(playa_alpha, playa_beta))

playa_sim_medians <- apply(playa_sims, 1, FUN = median)
playa_sim_means <- apply(playa_sims, 1, FUN = mean)
playa_sim_sds <- apply(playa_sims, 1, FUN = sd)

# Copper Mining Waste

copper_sims <- replicate(1e5, single_draw(copper_alpha, copper_beta))

copper_sim_medians <- apply(copper_sims, 1, FUN = median)
copper_sim_means <- apply(copper_sims, 1, FUN = mean)
copper_sim_sds <- apply(copper_sims, 1, FUN = sd)

# Paved road dust - Highway
highway_sims <- replicate(1e5, single_draw(highway_alpha, highway_beta))

highway_sim_medians <- apply(highway_sims, 1, FUN = median)
highway_sim_means <- apply(highway_sims, 1, FUN = mean)
highway_sim_sds <- apply(highway_sims, 1, FUN = sd)


# Excavation - Rock Crushing
rock_sims <- replicate(1e5, single_draw(rock_alpha, rock_beta))

rock_sim_medians <- apply(rock_sims, 1, FUN = median)
rock_sim_means <- apply(rock_sims, 1, FUN = mean)
rock_sim_sds <- apply(rock_sims, 1, FUN = sd)

# Cement
cement_sims <- replicate(1e5, single_draw(cement_alpha, cement_beta))

cement_sim_medians <- apply(cement_sims, 1, FUN = median)
cement_sim_means <- apply(cement_sims, 1, FUN = mean)
cement_sim_sds <- apply(cement_sims, 1, FUN = sd)


# Brake Wear
brake_sims <- replicate(1e5, single_draw(brake_alpha, brake_beta))

brake_sim_medians <- apply(brake_sims, 1, FUN = median)
brake_sim_means <- apply(brake_sims, 1, FUN = mean)
brake_sim_sds <- apply(brake_sims, 1, FUN = sd)


# Motor Vehicle Exhaust
exhaust_sims <- replicate(1e5, single_draw(exhaust_alpha, exhaust_beta))

exhaust_sim_medians <- apply(exhaust_sims, 1, FUN = median)
exhaust_sim_means <- apply(exhaust_sims, 1, FUN = mean)
exhaust_sim_sds <- apply(exhaust_sims, 1, FUN = sd)


# Baseline Rhizoplaca
baseline_sims <- replicate(1e5, single_draw(baseline_alpha, baseline_beta))

baseline_sim_medians <- apply(baseline_sims, 1, FUN = median)
baseline_sim_means <- apply(baseline_sims, 1, FUN = mean)
baseline_sim_sds <- apply(baseline_sims, 1, FUN = sd)


# Natural & anthropogenic
natural_sims <- replicate(1e5, single_draw(natural, empty_beta))

natural_sim_medians <- apply(natural_sims, 1, FUN = median)
natural_sim_means <- apply(natural_sims, 1, FUN = mean)
natural_sim_sds <- apply(natural_sims, 1, FUN = sd)

anthropogenic_sims <- replicate(1e5, single_draw(anthropogenic, empty_beta))

anthropogenic_sim_medians <- apply(anthropogenic_sims, 1, FUN = median)
anthropogenic_sim_means <- apply(anthropogenic_sims, 1, FUN = mean)
anthropogenic_sim_sds <- apply(anthropogenic_sims, 1, FUN = sd)


# Create matrices containing all simulated means, medians, and standard deviations for all profiles ##################################################################
######################################################################################################################################################################
sim_means_matrix <- rbind(baseline_sim_means, playa_sim_means, copper_sim_means, brake_sim_means, exhaust_sim_means, unpaved_dust_sim_means, coal_plant_sim_means)


sim_medians_matrix = rbind(baseline_sim_medians, playa_sim_medians, copper_sim_medians, 
                           brake_sim_medians, exhaust_sim_medians, unpaved_dust_sim_medians, 
                           coal_plant_sim_medians)

sim_sds_matrix = rbind(baseline_sim_sds, playa_sim_sds, copper_sim_sds, 
                           brake_sim_sds, exhaust_sim_sds, unpaved_dust_sim_sds, 
                           coal_plant_sim_sds)

# Match each column with its corresponding element
elem_names = c("Ca", "K",  "Mg", "N",  "P",  "S",  "Al", "As", "B", "Ba", "Cd", "Co", "Cr",
  "Cu", "Fe", "Mn", "Mo", "Na", "Ni", "Pb", "Se", "Si", "Sr", "Ti", "Zn")
colnames(sim_medians_matrix) = colnames(sim_sds_matrix) = colnames(sim_means_matrix) = elem_names


# Create table in LaTeX
#library("xtable")
#xtable(sim_medians_matrix[,order(elem_names)], digits=7)


