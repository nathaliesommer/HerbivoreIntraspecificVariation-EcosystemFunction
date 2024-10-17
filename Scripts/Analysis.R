### Scripts to run and objects needed ----

# Decomposition.R: combined_decomp
# SIR.R: SIR_final_data
# N-min.R: N_min_full_data
# Veg.R: functional_groups_wide, combined_diversity_long
# SoilPlantCN.R: CNdata


# Load necessary libraries
library(dplyr)
library(tidyr)

### Join data frames ----


# Prepare data from each script

# Decomposition data
combined_decomp <- combined_decomp %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Population, MassLoss, Year)

# SIR data
sir_data <- SIR_final_data %>%
  dplyr::select(Sample_ID, Year, CO2CperHourperg)

# N-mineralization data
n_min_data <- N_min_full_data %>%
  rename(Sample_ID = Sample.ID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, `Ammonium rate`, `Nitrate rate`, `Overall mineralization rate`)

# Vegetation biomass data
veg_biomass_data <- functional_groups_wide %>%
  rename(Sample_ID = Cage.ID) %>%
  dplyr::select(Sample_ID, Year, Population, Treatment, Transplant, Site, SORU_Biomass, POPRC_Biomass, MISC_Biomass)

# Plant diversity data
diversity_data <- combined_diversity_long %>%
  rename(Sample_ID = Cage.ID) %>%
  dplyr::select(Sample_ID, Year, PlantRichness, PlantDiversity)

# Soil and plant CN data
cn_data <- CNdata %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))

# Combine all data frames with inline type conversion
final_data <- combined_decomp %>%
  mutate(Year = as.character(Year)) %>%
  left_join(sir_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  left_join(n_min_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  left_join(veg_biomass_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  left_join(diversity_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  left_join(cn_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  # Select and rename the desired columns
  dplyr::select(
    Sample_ID,
    Year,
    Population = Population.x,  
    Site = Site.x,
    Treatment,
    Transplant,
    MassLoss,
    CO2CperHourperg,
    `Ammonium rate`,
    `Nitrate rate`,
    `Overall mineralization rate`,
    PlantRichness,
    PlantDiversity,
    SORU_Biomass,
    POPRC_Biomass,
    MISC_Biomass,
    PercentN_SOIL,
    PercentN_LITTER,
    PercentN_POPRC,
    PercentN_MISC,
    PercentN_SORU,
    PercentC_SOIL,
    PercentC_LITTER,
    PercentC_POPRC,
    PercentC_MISC,
    PercentC_SORU
  ) %>%
  mutate(
    Site_DailyMax = case_when(
      Site %in% c("YF", "DC", "SC") ~ "High",
      Site %in% c("UP", "FN") ~ "Low"
    )) %>% 
      mutate(Population_Origin = case_when(
        Population %in% c("YF", "DC", "SC") ~ "High",
        Population %in% c("UP", "FN") ~ "Low"
      )) %>% 
    mutate(Transplant_Temp = case_when(
      Transplant == "Home" ~ "Home",
      Transplant == "North" ~ "Low",
      Transplant == "South" ~ "High"
    ))


# Define potential response variables
response_vars <- c(
  "MassLoss", # y
  "CO2CperHourperg", # y
  "Ammonium rate",
  "Nitrate rate",
  "Overall mineralization rate", # y
  "SORU_Biomass", # y
  "POPRC_Biomass", # y
  "MISC_Biomass",
  "PercentN_SOIL", # y
  "PercentN_LITTER", # y
  "PercentN_POPRC", # y
  "PercentN_MISC",
  "PercentN_SORU", # y
  "PercentC_SOIL", # y
  "PercentC_LITTER",
  "PercentC_POPRC",
  "PercentC_MISC",
  "PercentC_SORU",
  "PlantRichness", # TO DO
  "PlantDiversity" # TO DO
)

# Create a new data object that calculates the difference between 2021 and 2023 for each response variable
final_data_diff <- final_data %>%
  filter(Year %in% c(2021, 2023)) %>%
  group_by(Sample_ID, Population, Site, Treatment, Transplant, Year) %>%
  summarise(across(all_of(response_vars), mean, .names = "{col}"), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = all_of(response_vars), names_sep = "_") %>%
  mutate(
    Diff_MassLoss = MassLoss_2023 - MassLoss_2021,
    Diff_CO2CperHourperg = CO2CperHourperg_2023 - CO2CperHourperg_2021,
    Diff_Ammonium_rate = `Ammonium rate_2023` - `Ammonium rate_2021`,
    Diff_Nitrate_rate = `Nitrate rate_2023` - `Nitrate rate_2021`,
    Diff_Overall_mineralization_rate = `Overall mineralization rate_2023` - `Overall mineralization rate_2021`,
    Diff_SORU_Biomass = SORU_Biomass_2023 - SORU_Biomass_2021,
    Diff_POPRC_Biomass = POPRC_Biomass_2023 - POPRC_Biomass_2021,
    Diff_MISC_Biomass = MISC_Biomass_2023 - MISC_Biomass_2021,
    Diff_PercentN_SOIL = PercentN_SOIL_2023 - PercentN_SOIL_2021,
    Diff_PercentN_LITTER = PercentN_LITTER_2023 - PercentN_LITTER_2021,
    Diff_PercentN_POPRC = PercentN_POPRC_2023 - PercentN_POPRC_2021,
    Diff_PercentN_MISC = PercentN_MISC_2023 - PercentN_MISC_2021,
    Diff_PercentN_SORU = PercentN_SORU_2023 - PercentN_SORU_2021,
    Diff_PercentC_SOIL = PercentC_SOIL_2023 - PercentC_SOIL_2021,
    Diff_PercentC_LITTER = PercentC_LITTER_2023 - PercentC_LITTER_2021,
    Diff_PercentC_POPRC = PercentC_POPRC_2023 - PercentC_POPRC_2021,
    Diff_PercentC_MISC = PercentC_MISC_2023 - PercentC_MISC_2021,
    Diff_PercentC_SORU = PercentC_SORU_2023 - PercentC_SORU_2021,
    Diff_PlantRichness = PlantRichness_2023 - PlantRichness_2021,
    Diff_PlantDiversity = PlantDiversity_2023 - PlantDiversity_2021,
    Population_Origin = case_when(
      Population %in% c("YF", "DC", "SC") ~ "High",
      Population %in% c("UP", "FN") ~ "Low"
    ),
    Transplant_Temp = case_when(
      Transplant == "Home" ~ "Home",
      Transplant == "North" ~ "Low",
      Transplant == "South" ~ "High",
    )
  ) %>%
  dplyr::select(Sample_ID, Population, Site, Treatment, Transplant, Population_Origin, Transplant_Temp, starts_with("Diff_"))

# Relevel for easier interpretation
final_data_diff$Treatment <- as.factor(final_data_diff$Treatment)
final_data_diff$Treatment <- relevel(final_data_diff$Treatment, ref = "Vegetation")
final_data_diff$Transplant_Temp <- as.factor(final_data_diff$Transplant_Temp)
final_data_diff$Transplant_Temp <- relevel(final_data_diff$Transplant_Temp, ref = "Home")



# Model Structure: Difference in Response ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site)


### Models ----

library(lme4)
library(ggplot2)
library(performance)
library(DHARMa)
library(car) # for VIF calculation
library(boot)
library(MuMIn)

#### Decomposition ----

hist(final_data_diff$Diff_MassLoss)
Decomp_model <- lmer(Diff_MassLoss ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site), 
                      data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(Decomp_model)
simulation_output <- simulateResiduals(fittedModel = Decomp_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)
# significant KS but no overdispersion

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_MassLoss ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values) 


# LMM for Decomp_model is the final model. LMMs are somewhat robust to deviations from normality. 
# Instead, use non-parametric bootstrapping to generate confidence intervals for the fixed effects 
                                            # (does not rely on normality assumption in residuals)

summary(Decomp_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(Decomp_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(Decomp_model))
fixed_effects <- fixef(Decomp_model)
decomp_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(decomp_effect_summary)

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(Decomp_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(Decomp_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(Decomp_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])















#### SIR ----

hist(final_data_diff$Diff_CO2CperHourperg)

SIR_model <- lmer(Diff_CO2CperHourperg ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(SIR_model)
simulation_output <- simulateResiduals(fittedModel = SIR_model, n = 1000)
plot(simulation_output) # heteroscedasticity in the model
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_CO2CperHourperg ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)

# LMM for SIR_model is the final model. The heteroscedasticity is not severe, and there is no indication of overdispersion
# Use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(SIR_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(SIR_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(SIR_model))
fixed_effects <- fixef(SIR_model)
SIR_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(SIR_effect_summary)

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(SIR_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(SIR_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(SIR_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])



















#### SORU Biomass ----

hist(final_data_diff$Diff_SORU_Biomass)

SORU_biomass_model <- lmer(Diff_SORU_Biomass ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(SORU_biomass_model)
plot(simulateResiduals(fittedModel = SORU_biomass_model, n = 1000)) #significant KS, mild heteroscedasticity
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_SORU_Biomass ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)


# Use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(SORU_biomass_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(SORU_biomass_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(SORU_biomass_model))
fixed_effects <- fixef(SORU_biomass_model)
SORUbiomass_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(SORUbiomass_effect_summary) # interpret

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(SORU_biomass_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(SORU_biomass_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(SORU_biomass_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])










#### POPRC_Biomass ----
hist(final_data_diff$Diff_POPRC_Biomass)

POPRC_biomass_model <- lmer(Diff_POPRC_Biomass ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(POPRC_biomass_model)
simulation_output <- simulateResiduals(fittedModel = POPRC_biomass_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_POPRC_Biomass ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)



# Diagnostics look great; still use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(POPRC_biomass_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(POPRC_biomass_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(POPRC_biomass_model))
fixed_effects <- fixef(POPRC_biomass_model)
POPRCbiomass_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(POPRCbiomass_effect_summary) # interpret

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(POPRC_biomass_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(POPRC_biomass_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(POPRC_biomass_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])












#### N-min rate ----

Nmin_model <- lmer(Diff_Overall_mineralization_rate ~ Treatment + Population_Origin * Transplant_Temp + (1 | Site), data = final_data_diff)
plot(Nmin_model)

# 1. Diagnostics - DHARMa package
hist(final_data_diff$Diff_Overall_mineralization_rate)
simulation_output <- simulateResiduals(fittedModel = Nmin_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) 

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_Overall_mineralization_rate ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)



# Diagnostics look good, but we will still use non-parametric bootstrapping 

summary(Nmin_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(Nmin_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(Nmin_model))
fixed_effects <- fixef(Nmin_model)
Nmin_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(Nmin_effect_summary) # interpret
# two significant parameters, but the effect size is negligible

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(Nmin_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(Nmin_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(Nmin_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])











#### TO DO: Soil %C ----
SoilC_model <- lmer(Diff_PercentC_SOIL ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)
plot(SoilC_model)

# 1. Diagnostics - DHARMa package
hist(final_data_diff$Diff_PercentC_SOIL)
simulation_output <- simulateResiduals(fittedModel = SoilC_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through ChatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PercentC_SOIL ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)









#### TO DO: Soil %N ----
hist(final_data_diff$Diff_PercentN_SOIL)
SoilN_model <- lmer(Diff_PercentN_SOIL ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(SoilN_model)
simulation_output <- simulateResiduals(fittedModel = SoilN_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through ChatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PercentN_SOIL ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)












#### TO DO: Litter %N ----
hist(final_data_diff$Diff_PercentN_LITTER)
LitterN_model <- lmer(Diff_PercentN_LITTER ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(LitterN_model)
simulation_output <- simulateResiduals(fittedModel = LitterN_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through ChatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PercentN_LITTER ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)






#### TO DO: POPRC %N ----

hist(final_data_diff$Diff_PercentN_POPRC)
POPRCN_model <- lmer(Diff_PercentN_POPRC ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(POPRCN_model)
simulation_output <- simulateResiduals(fittedModel = POPRCN_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through ChatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PercentN_POPRC ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)






#### TO DO: SORU %N ----
hist(final_data_diff$Diff_PercentN_SORU)
SORUN_model <- lmer(Diff_PercentN_SORU ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(SORUN_model)
simulation_output <- simulateResiduals(fittedModel = SORUN_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through ChatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PercentN_SORU ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)











#### TO DO: PlantRichness ----
hist(final_data_diff$Diff_PlantRichness)
PlantRichness_model <- lmer(Diff_PlantRichness ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(PlantRichness_model)
simulation_output <- simulateResiduals(fittedModel = PlantRichness_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # zero-inflated, run through chatGPT

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PlantRichness ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)








#### PlantDiversity ----
hist(final_data_diff$Diff_PlantDiversity)
PlantDiversity_model <- lmer(Diff_PlantDiversity ~ Treatment + Population_Origin * Transplant_Temp + (1| Site), data = final_data_diff)

# 1. Diagnostics - DHARMa package
plot(PlantDiversity_model)
simulation_output <- simulateResiduals(fittedModel = PlantDiversity_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(Diff_PlantDiversity ~ Treatment + Population_Origin * Transplant_Temp, data = final_data_diff)
vif_values <- vif(fixed_model)
print(vif_values)

# Diagnostics look good, but we will still use non-parametric bootstrapping 

summary(PlantDiversity_model)

##### Fixed effects and confidence intervals together ----
boot_model <- bootMer(PlantDiversity_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(PlantDiversity_model))
fixed_effects <- fixef(PlantDiversity_model)
PlantDiversity_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(PlantDiversity_effect_summary) # interpret
# two significant parameters, but the effect size is negligible

##### Variance partitioning ----

# Extract the variance components of the random effects
random_effects_variance <- as.data.frame(VarCorr(PlantDiversity_model))

# Extract the variance for the random effects and residuals
site_variance <- random_effects_variance$vcov[1]  
residual_variance <- attr(VarCorr(PlantDiversity_model), "sc")^2

# Calculate total variance
total_variance <- site_variance + residual_variance

# Calculate the proportion of variance explained by the random effect (Site)
random_effect_variance_proportion <- site_variance / total_variance

# Calculate the proportion of variance explained by the residuals
residual_variance_proportion <- residual_variance / total_variance

# Print results
cat("Proportion of variance explained by Site (Random Effect):", random_effect_variance_proportion)
cat("Proportion of variance explained by Residuals:", residual_variance_proportion)

# Calculate the marginal and conditional R² values
r_squared <- r.squaredGLMM(Nmin_model)

# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])




