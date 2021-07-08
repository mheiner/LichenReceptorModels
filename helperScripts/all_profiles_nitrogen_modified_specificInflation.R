# ALL PROFILES, nitrogen modified based on Nitrate/Ammonium concentrations

# Exclude the 5 elements we're not using: (Br, Cl, F, Rb, V)
# indices of these elements:
removed_elements <- c(11,12,17,24,29)
element_names <- c("Ca", "K",  "Mg", "N",  "P",  "S",  "Al", "As", "B", "Ba", "Cd", "Co", "Cr",
  "Cu", "Fe", "Mn", "Mo", "Na", "Ni", "Pb", "Se", "Si", "Sr", "Ti", "Zn")

#####################################################################################################################################
# Contructing the alpha and beta prior matrices ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Here, we use "alpha" and "beta", which are used in independent gamma priors which, when normalized,
# create draws from the "Generalized Dirichlet" prior that we are using on the source profiles.

# In this case, alpha (shape) = mean^2 / sd^2, and beta (rate) = mean / sd^2

# Mean and standard deviation values were obtained from SPECIATE, existing literature, or previous estimation (Baseline)
# Each profile will follow the general pattern outlined below:

# Create a vector of means
# Create a vector of standard deviations
# Adjust the standard deviation by some factor
  # In the case where all standard deviations are missing, insert the mean value for the same element,
  # adjusted by some "A_factor", then continue.
# Remove the elements we are not considering in the analysis
# Replace missing mean values (NA's) with the raw mean of the other mean values
# When a mean value should be 0, replace with 1e-6 to facilitate sampling later on
# For standard deviations with missing values, use the raw mean multiplied by some "A_factor"
# Compute alpha and beta based on the calculation above


#####################################################################################################################################

# Use inflation factors to adjust scaling of standard deviations
inflation_factor <- 1
A_factor <- 1

# Agriculture Soil: Nitrate = 0, sd = 0.6733. Ammonium = 0, sd = 0.1872.
# I'll input "0" for the mean for nitrogen, and put the sd as 0.5 (keeping sd
# at NA results in a sd of 1.922 which seems far too high)
agriculture_means <- c(1.678,	2.488,	NA,	0,	0.1518,	0.0633,	7.4502,	0,	NA,	0,	0,	0.118,	0,	0,	0,	0.0371,	NA,	4.7058,	0.1186,	0,
                       NA,	0,	0,	0.0163,	0,	25.0138,	0.0235,	0.5193,	0,	0.0464)
agriculture_sds <- c(0.3177,	0.3959,	NA,	0.5,	0.0382,	0.0434,	1.139,	0.0283,	NA,	0.7295,	0.0118,	0.0896,	0.1179,	0.0705,	0.0228,	0.027,
                     NA,	0.6402,	0.0281,	0.0436,	NA,	0.0086,	0.0427,	0.0142,	0.0143,	4.337,	0.0167,	0.1229,	0.0981,	0.0203)
agriculture_sds <- agriculture_sds * inflation_factor


agriculture_means <- agriculture_means[-removed_elements]
agriculture_sds <- agriculture_sds[-removed_elements]

agriculture_M <- mean(agriculture_means, na.rm=TRUE)

agriculture_means[is.na(agriculture_means)] <- agriculture_M
agriculture_means[agriculture_means == 0] <- 1e-6
agriculture_sds[is.na(agriculture_sds)] <- agriculture_M * A_factor

agriculture_alpha <- agriculture_means^2 / agriculture_sds^2
agriculture_beta <- agriculture_means / agriculture_sds^2

# Coal Combustion: No ammonium. Nitrate = 0, sd = 0.2. So, I'll put a mean on
# nitrogen at 0 with a sd of 0.2 to be safe.
coal_combustion_means <- c(3.4536,	0.4644,	NA,	0,	0.9372,	2.948,	5.968,	0,	NA,	1.3315,	0.0147,	0.0629,	0,	NA,	0.0176,	0.0179,	NA,	2.916,
                           0.0284,	0,	NA,	0.0072,	0.068,	0.0053,	0.0406,	9.0112,	0.1964,	0.4315,	0,	0.0797)
coal_combustion_sds <- c(1.0411,	0.0602,	NA,	0.2,	0.6322,	2.729,	0.5247,	0.0164,	NA,	1.0801,	0.0154,	0.0221,	0.0341,	NA,	0.0041,	0.0112,	NA,
                         0.3827,	0.0139,	0.0134,	NA,	0.0019,	0.0336,	0.0043,	0.0407,	0.5675,	0.0686, 0.0651,	0.0734,	0.0341)
c_c_means <- coal_combustion_means
c_c_sds <- coal_combustion_sds
c_c_sds <- c_c_sds * inflation_factor

c_c_means <- c_c_means[-removed_elements]
c_c_sds <- c_c_sds[-removed_elements]

c_c_M <- mean(coal_combustion_means, na.rm=TRUE)

c_c_means[is.na(c_c_means)] <- c_c_M
c_c_means[c_c_means == 0] <- 1e-6
c_c_sds[is.na(c_c_sds)] <- c_c_M * A_factor

c_c_alpha <- c_c_means^2 / c_c_sds^2
c_c_beta <- c_c_means / c_c_sds^2

# Unpaved Road Dust: nitrate = 0, sd=0.6731. Ammonium = 0, sd = 0.1317
# I'll use a mean on Nitrogen of 0 with a sd of 0.5 (similar to above in Agriculture)
unpaved_means <- c(2.163,	2.8299,	NA,	0,	0.1603,	0.2808,	7.4822,	0,	NA,	0,	0,	0.1519,	0,	0,	0.0313,	0.0474,	NA,	5.5128,	0.1372,	0,	NA,
                   0.0091,	0.0288,	0.0184,	0,	24.2969,	0.0312,	0.5258,	0,	0.0502)
unpaved_sds <- c(1.0444,	0.4949,	NA,	0.5,	0.044,	0.3884,	0.9315,	0.0226,	NA,	0.5473,	0.0078,	0.0755,	0.0881,	0.0869,	0.0161,	0.0307,	NA,	2.1152,
                 0.0509,	0.0331,	NA,	0.0057,	0.0284,	0.0093,	0.0108,	4.0089,	0.0112,	0.1289,	0.0646,	0.021)

unpaved_adjust = 2.0
unpaved_sds <- unpaved_sds * inflation_factor * unpaved_adjust

unpaved_means <- unpaved_means[-removed_elements]
unpaved_sds <- unpaved_sds[-removed_elements]

unpaved_M <- mean(unpaved_means, na.rm=TRUE)

unpaved_means[is.na(unpaved_means)] <- unpaved_M
unpaved_means[unpaved_means == 0] <- 1e-6
unpaved_sds[is.na(unpaved_sds)] <- unpaved_M * A_factor

unpaved_alpha <- unpaved_means^2 / unpaved_sds^2
unpaved_beta <- unpaved_means / unpaved_sds^2

# Coal-Fired Power Plant: No nitrate/nitrite/ammonium listed, so keep Nitrogen as NA
coal_plant_means <- c(1.393,	1.379,	NA,	NA,	0.427,	1.803,	14.843,	0.055,	NA,	0.056,	NA,	0.099,	0.005,	NA,	0.054,	0.029,	NA,	9.077,
                      0.045,	NA,	NA,	0.039,	0.03,	0.011,	0.018,	23.29,	0.144,	0.99,	0.072,	0.055)
coal_plant_sds <- c(0.129,	0.127,	NA,	NA,	0.07,	0.221,	1.377,	0.006,	NA,	0.046,	NA,	0.02,	0.006,	NA,	0.005,	0.003,	NA,	0.831,	0.004,
                    NA,	NA,	0.004,	0.004,	0.001,	0.003,	2.15,	0.013,	0.091,	0.007,	0.006)
coal_plant_sds <- coal_plant_sds * inflation_factor

coal_plant_means <- coal_plant_means[-removed_elements]
coal_plant_sds <- coal_plant_sds[-removed_elements]

coal_plant_M <- mean(coal_plant_means, na.rm=TRUE)

coal_plant_means[is.na(coal_plant_means)] <- coal_plant_M
coal_plant_means[coal_plant_means == 0] <- 1e-6
coal_plant_sds[is.na(coal_plant_sds)] <- coal_plant_M * A_factor

coal_plant_alpha <- coal_plant_means^2 / coal_plant_sds^2
coal_plant_beta <- coal_plant_means / coal_plant_sds^2

# Oil Refinery: Nitrate = 0.1089, sd = 0.1234.
# Try nitrogen = (14/62)*0.1089 = 0.0246, with sd = 0.1
oil_means <- c(2.1993,	1.7413,	0.2768,	0.0246,	0.074,	0.1318,	9.4368,	0.0022,	NA,	0.0662,	0.0008,	0.0632,	0,	0.002,	0.0332,	0.0127,	NA,
               2.8723,	0.0564,	0.006,	0.4126,	0.0939,	0.0218,	0.0092,	0,	26.5176,	0.019,	0.388,	0.1685,	0.0906)
oil_sds <- c(0.4095,	0.2993,	0.2265, 0.1,	0.0306,	0.0363,	1.6607,	0.0123,	NA,	0.3367,	0.0057,	0.0478,	0.0568,	0.052,	0.0114,	0.006,	NA,
             0.2493,	0.0083,	0.0144,	0.1954,	0.0114,	0.0162,	0.0051,	0.006,	4.8632,	0.0045,	0.101,	0.0454,	0.0089)
oil_sds <- oil_sds * inflation_factor

oil_means <- oil_means[-removed_elements]
oil_sds <- oil_sds[-removed_elements]

oil_M <- mean(oil_means, na.rm=TRUE)

oil_means[is.na(oil_means)] <- oil_M
oil_means[oil_means == 0] <- 1e-6
oil_sds[is.na(oil_sds)] <- oil_M * A_factor

oil_alpha <- oil_means^2 / oil_sds^2
oil_beta <- oil_means / oil_sds^2

# Avg. Fine Playa:
playa_means <- c(124646,	4067,	24433,	1000,	5000,	20000,	2620,	14.2,	148.5,	144,	NA,	NA,	0.24,	1e-2,	19.6,	20.3,	NA,	3936,	153,	20.39,	96378,
                 1e-2, 10,	15.24,	0.34,	30000,	1320,	89.7,	8.97,	1)
playa_sds <- c(73976,	4392,	10013,	5000,	4000,	10000,	1619,	9.2,	140.8,	76,	NA,	NA,	0.16,	1e-1,	75.2,	22,	NA,	2919,	93,	54.66,	108002,	1e-1,	8,
               10.25,	0.46,	10000,	945,	64.3,	4.73,	.5)

playa_adjust = 0.5
playa_sds <- playa_sds * inflation_factor * playa_adjust

playa_means <- playa_means[-removed_elements]
playa_sds <- playa_sds[-removed_elements]

playa_M <- mean(playa_means, na.rm=TRUE)

playa_means[is.na(playa_means)] <- playa_M
playa_sds[is.na(playa_sds)] <- playa_M * A_factor

playa_alpha <- playa_means^2 / playa_sds^2
playa_beta <- playa_means / playa_sds^2

# Copper Mining Waste: no nitrate/nitrite/ammonium, keep Nitrogen as NA
copper_means <- c(1.687,	2.155,	NA,	NA,	0.382,	3.354,	10.263,	0.045,	NA,	NA,	0.005,	0.096,	0.009,	NA,	0.041,	0.296,	NA,	10.658,
                  0.102,	NA,	NA,	NA,	0.799,	0.01,	NA,	24.707,	NA,	0.485,	0.038,	0.682)
copper_adjust = 1.0

copper_A <- 0.5
copper_sds <- copper_A * copper_means

copper_sds <- copper_sds * inflation_factor * copper_adjust

copper_means <- copper_means[-removed_elements]
copper_sds <- copper_sds[-removed_elements]

copper_M <- mean(copper_means, na.rm=TRUE)

copper_means[is.na(copper_means)] <- copper_M
copper_sds[is.na(copper_sds)] <- copper_M * A_factor

copper_alpha <- copper_means^2 / copper_sds^2
copper_beta <- copper_means / copper_sds^2


# Paved Road Dust - Highway: Nitrate=0.0699, sd = 0.0328. Ammonium=0.0776, sd=0.0329.
# Use: Nitrogen mean = (14/62)*0.0699 + (14/18)*0.0776 = 0.076. Try sd = 0.1
highway_means <- c(4.9452,	2.2023,	0.7336,	0.076,	0.0631,	1.2905,	7.9698,	0,	NA,	0.1081,	0.0012,	0.4855,	0,
                   0.0007,	0.0068,	0.0048,	NA,	3.1407,	0.0551,	0.0007,	0.4415,	0.0087,	0.014,	0.0097,
                   0,	24.3753,	0.043,	0.3073,	0.0169,	0.0418)
highway_sds <- c(0.868,	0.4519,	0.0636,	0.1,	0.0299,	0.0921,	2.4017,	0.0041,	NA,	0.0312,	0.0019,	0.1447,
                 0.0169,	0.0476,	0.002,	0.0008,	NA,	0.2229,	0.0046,	0.0041,	0.0633,	0.0008,	0.0021,
                 0.0009,	0.0017,	7.8115,	0.0032,	0.0253,	0.0198,	0.003)
highway_sds <- highway_sds * inflation_factor

highway_means <- highway_means[-removed_elements]
highway_sds <- highway_sds[-removed_elements]

highway_M <- mean(highway_means, na.rm=TRUE)

highway_means[is.na(highway_means)] <- highway_M
highway_means[highway_means == 0] <- 1e-6
highway_sds[is.na(highway_sds)] <- highway_M * A_factor

highway_alpha <- highway_means^2 / highway_sds^2
highway_beta <- highway_means / highway_sds^2

# Excavation - Rock Crushing: no nitrate/nitrite/ammonium, keep Nitrogen = NA
rock_means <- c(1.732,	2.721,	NA,	NA,	NA,	0,	9.681,	0,	NA,	0.027,	NA,	0,	NA,	NA,	0.007,	0.005,
                NA,	2.908,	0.117,	NA,	NA,	0.003,	0.014,	0.012,	NA,	29.414,	0.057,	0.339,	0.008, 0.01)
rock_sds <- c(0.087,	0.136,	NA,	NA,	NA,	NA,	1.452,	0.001,	NA,	0.027,	NA,	NA,	NA,	NA,	0.001,	0.001,
              NA,	0.147,	0.006,	NA,	NA,	0.001,	0.002,	0.001,	NA,	4.412,	0.003,	0.022,	0.006, 0.01)
rock_sds <- rock_sds * inflation_factor

rock_means <- rock_means[-removed_elements]
rock_sds <- rock_sds[-removed_elements]

rock_M <- mean(rock_means, na.rm=TRUE)

rock_means[is.na(rock_means)] <- rock_M
rock_means[rock_means == 0] <- 1e-6
rock_sds[is.na(rock_sds)] <- rock_M * A_factor

rock_alpha <- rock_means^2 / rock_sds^2
rock_beta <- rock_means / rock_sds^2


# Cement: nitrate=0, sd=0.0141. Ammonium=0.0033, sd=0.0141
# Use Nitrogen mean = 0 + (14/18)*0.0033 = 0.0026. Try sd = 0.01
cement_means <- c(36.8033,	0.7201,	0.5156,	0.0026,	0.1633,	2.1793,	3.8864,	0.0004,	NA,	0.179,	0,	0,	0,
                  0,	0.0094,	0.0049,	NA,	2.398,	0.0273,	0,	0,	0.0095,	0.0061,	0.0038,	0.0008,	11.8661,
                  0.1255,	0.1379,	0.0151,	0.0128)
cement_sds <- c(6.4548,	0.1712,	0.0467,	0.01,	0.0758,	0.1543,	1.1722,	0.002,	NA,	0.0205,	0.001,	0.2355,	0.011,
                0.0363,	0.0015,	0.0004,	NA,	0.1697,	0.0024,	0.0019,	0.1149,	0.0014,	0.0009,	0.0004,	0.0008,
                3.8022,	0.0089,	0.0134,	0.0045,	0.001)
cement_sds <- cement_sds * inflation_factor

cement_means <- cement_means[-removed_elements]
cement_sds <- cement_sds[-removed_elements]

cement_M <- mean(cement_means, na.rm=TRUE)

cement_means[is.na(cement_means)] <- cement_M
cement_means[cement_means == 0] <- 1e-6
cement_sds[is.na(cement_sds)] <- cement_M * A_factor

cement_alpha <- cement_means^2 / cement_sds^2
cement_beta <- cement_means / cement_sds^2


# Brake Wear: ammonium = 0, Nitrate = 0.16
# Try Nitrogen mean = (14/62)*0.16 = 0.036. Leave sd=NA like all the others for simplicity
brake_means <- c(0.11,	0.02,	8.3,	0.0334,	0,	1.28,	0.03,	0,	NA,	5.45,	NA,	0.15,	NA,	NA,	0.12,	1.15,	NA,
                 28.7,	0.17,	0.37,	0,	0.07,	0.01,	0.01,	0,	6.79,	0.07,	0.36,	0.07,	0.03)
brake_A <- 0.5
brake_means[ brake_means == 0] <- 1e-6

brake_sds <- brake_A * brake_means

brake_sds <- brake_sds * inflation_factor

brake_means <- brake_means[-removed_elements]
brake_sds <- brake_sds[-removed_elements]

brake_M <- mean(brake_means, na.rm=TRUE)

brake_means[is.na(brake_means)] <- brake_M

brake_sds[is.na(brake_sds)] <- brake_M * A_factor

brake_alpha <- brake_means^2 / brake_sds^2
brake_beta <- brake_means / brake_sds^2

# Motor Vehicle Exhaust: nitrate=0.06, sd=NA
# Try Nitrogen mean = (14/62)*(0.06) = 0.0135. Keep sd=NA.
exhaust_means <- c(0.42,	0.01,	0.14,	NA,	0.27,	0.37,	0.1,	NA,	NA,	0.01,	NA,	0.39,	NA,	NA,	0.01,	0.02,
                   NA,	1.27,	0.01,	NA,	0.01,	0.01,	0.08,	NA,	NA,	1.61,	NA,	NA,	NA,	0.49)
exhaust_A <- 0.5
exhaust_sds <- exhaust_A * exhaust_means
exhaust_sds <- exhaust_sds * inflation_factor

exhaust_means <- exhaust_means[-removed_elements]
exhaust_sds <- exhaust_sds[-removed_elements]

exhaust_M <- mean(exhaust_means, na.rm=TRUE)

exhaust_means[is.na(exhaust_means)] <- exhaust_M
exhaust_sds[is.na(exhaust_sds)] <- exhaust_M * A_factor

exhaust_alpha <- exhaust_means^2 / exhaust_sds^2
exhaust_beta <- exhaust_means / exhaust_sds^2

# Baseline Rhizoplaca
baseline_rhizoplaca_means <- c(0.577469279,	0.081169754,	0.018432954,	0.220858918,	0.023583944,
                               0.024024262,	0.007371760,	0.000042000,	0.000157807,	0.000516918,	NA,
                               NA,	0.000017700,	0.000019100,	0.000083200,	0.000125271,	NA,
                               0.037699008,	0.000846568,	0.000003560,	0.002550509,	0.000033600,
                               0.000209994,	NA,	0.000002880,	0.001839548,	0.000943357,	0.001376008,
                               NA,	0.000622160)
baseline_rhizoplaca_sds <- c(0.057221061,	0.035932365,	0.006025406,	0.066461972,	0.011731145,
                             0.009622388,	0.028655619,	0.000026100,	0.000113140,	0.000229714,
                             NA,	NA,	0.000008610,	0.000012000,	0.000031700,	0.000040000,	NA,
                             0.012846014,	0.000356402,	0.000009690,	0.001376478,	0.000019000,
                             0.000059200,	NA,	0.000009960,	0.001849106,	0.000236535,	0.000656630,
                             NA,	0.000222225)
baseline_adjust = 3.0
baseline_rhizoplaca_sds <- baseline_rhizoplaca_sds * inflation_factor * baseline_adjust

baseline_rhizoplaca_means <- baseline_rhizoplaca_means[-removed_elements]
baseline_rhizoplaca_sds <- baseline_rhizoplaca_sds[-removed_elements]

baseline_M <- mean(baseline_rhizoplaca_means, na.rm=TRUE)

baseline_rhizoplaca_means[is.na(baseline_rhizoplaca_means)] <- baseline_M
baseline_rhizoplaca_sds[is.na(baseline_rhizoplaca_sds)] <- baseline_M * A_factor

baseline_alpha <- baseline_rhizoplaca_means^2 / baseline_rhizoplaca_sds^2
baseline_beta <- baseline_rhizoplaca_means / baseline_rhizoplaca_sds^2

# Assemble alpha and beta matrices containing alpha and beta values for each profile

alpha_matrix <- rbind(baseline_alpha, playa_alpha, copper_alpha, brake_alpha, exhaust_alpha, unpaved_alpha, coal_plant_alpha, cement_alpha)
beta_matrix <- rbind(baseline_beta, playa_beta, copper_beta, brake_beta, exhaust_beta, unpaved_beta, coal_plant_beta, cement_beta)

# Label each column with its associated element
colnames(alpha_matrix) = colnames(beta_matrix) = element_names

alpha_matrix

# Create weak "natural" and "anthropogenic" source profile priors

empty_alpha <- rep(1/ncol(alpha_matrix), ncol(alpha_matrix))

natural_weights <- c(5, 5, 5, 1, 1, 1, 1, 1, 1, 1, 0.2, 0.2, 0.2, 0.2, 1, 1, 0.2, 5, 0.2, 0.2, 0.2,
                     1, 5, 1, 0.2)
anthropogenic_weights <- c(0.2, 0.2, 0.2, 0.2, 1, 1, 5, 1, 1, 1, 5, 5, 5, 5, 1, 1, 5, 0.2, 5, 5, 
                           5, 5, 1, 1, 5)

natural <- empty_alpha * natural_weights * 4
anthropogenic <- empty_alpha * anthropogenic_weights * 4
empty_beta <- rep(1, ncol(alpha_matrix))


