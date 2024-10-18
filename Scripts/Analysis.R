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
decomp_data <- combined_decomp %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Population, MassLoss, Year)

sir_data <- SIR_final_data %>%
  mutate(Population = substr(Sample_ID, 1, 2)) %>%  # Add Population column
  dplyr::select(Sample_ID, Year, Population, CO2CperHourperg)

n_min_data <- N_min_full_data %>%
  rename(Sample_ID = Sample.ID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, `Ammonium rate`, `Nitrate rate`, `Overall mineralization rate`)

veg_biomass_data <- functional_groups_wide %>%
  rename(Sample_ID = Cage.ID) %>%
  dplyr::select(Sample_ID, Year, Population, Treatment, Transplant, Site, SORU_Biomass, POPRC_Biomass, MISC_Biomass)

diversity_data <- combined_diversity_long %>%
  rename(Sample_ID = Cage.ID) %>%
  mutate(Population = substr(Sample_ID, 1, 2)) %>%  # Add Population column
  dplyr::select(Sample_ID, Year, Population, PlantRichness, PlantDiversity)

cn_data <- CNdata %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))

# Define potential response variables
response_vars <- c(
  "MassLoss",
  "CO2CperHourperg",
  "Ammonium rate",
  "Nitrate rate",
  "Overall mineralization rate",
  "SORU_Biomass",
  "POPRC_Biomass",
  "MISC_Biomass",
  "PercentN_SOIL",
  "PercentN_LITTER",
  "PercentN_POPRC",
  "PercentN_MISC",
  "PercentN_SORU",
  "PercentC_SOIL",
  "PercentC_LITTER",
  "PercentC_POPRC",
  "PercentC_MISC",
  "PercentC_SORU",
  "PlantRichness",
  "PlantDiversity"
)

# Combine all data frames
combined_data <- decomp_data %>%
  mutate(Year = as.character(Year)) %>%
  left_join(sir_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(n_min_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(veg_biomass_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(diversity_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(cn_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  # Ensure Site column is present
  mutate(Site = coalesce(Site.x, Site.y, Site)) %>%
  # Remove duplicates
  distinct(Sample_ID, Year, .keep_all = TRUE)

# Select and mutate final data
final_data <- combined_data %>%
  dplyr::select(
    Sample_ID,
    Year,
    Population,  
    Site,
    Treatment,
    Transplant,
    all_of(response_vars)
  ) %>%
  mutate(
    Population_Origin = case_when(
      Population %in% c("YF", "DC", "SC") ~ "High",
      Population %in% c("UP", "FN") ~ "Low"
    ),
    Transplant_Temp = case_when(
      Transplant == "Home" ~ "Home",
      Transplant == "North" ~ "Low",
      Transplant == "South" ~ "High"
    )
  )

# Create final_data_year with variables spread by Year
final_data_year <- final_data %>%
  pivot_wider(names_from = Year, values_from = all_of(response_vars), names_sep = "_") %>%
  mutate(
    Site_DailyMax = case_when(
      Site %in% c("YF", "DC", "SC") ~ "High",
      Site %in% c("UP", "FN") ~ "Low"
    )
  )

# Create final_data_diff for differences between years
final_data_diff <- final_data %>%
  filter(Year %in% c("2021", "2023")) %>%
  pivot_wider(names_from = Year, values_from = all_of(response_vars), names_sep = "_") %>%
  mutate(across(starts_with("2023_"), 
                ~ . - get(sub("2023_", "2021_", cur_column())), 
                .names = "Diff_{.col}")) %>%
  dplyr::select(Sample_ID, Population, Site, Treatment, Transplant, Population_Origin, Transplant_Temp, starts_with("Diff_"))


# Relevel for easier interpretation
final_data_diff$Treatment <- as.factor(final_data_diff$Treatment)
final_data_diff$Treatment <- relevel(final_data_diff$Treatment, ref = "Vegetation")
final_data_diff$Transplant_Temp <- as.factor(final_data_diff$Transplant_Temp)
final_data_diff$Transplant_Temp <- relevel(final_data_diff$Transplant_Temp, ref = "Home")

final_data_year$Treatment <- as.factor(final_data_year$Treatment)
final_data_year$Treatment <- relevel(final_data_year$Treatment, ref = "Vegetation")
final_data_year$Transplant <- as.factor(final_data_year$Transplant)
final_data_year$Transplant <- relevel(final_data_year$Transplant, ref = "Home")

# Model Structure: 2023 Response ~ 2021 Baseline + Treatment + Population_Origin * Transplant_Temp + (1 | Site)


### LMMs ----

library(lme4)
library(ggplot2)
library(performance)
library(DHARMa)
library(car) # for VIF calculation
library(boot)
library(MuMIn)

#### Decomposition ----

hist(final_data_year$MassLoss_2023)
Decomp_model <- lmer(MassLoss_2023 ~ MassLoss_2021 + Treatment + Population_Origin * Transplant_Temp + 
                       (1 | Site), 
                      data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(Decomp_model)
simulation_output <- simulateResiduals(fittedModel = Decomp_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)
# significant KS but no overdispersion

# 2. Multicollinearity Assessment
fixed_model <- lm(MassLoss_2023 ~ MassLoss_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 


# LMM for Decomp_model is the final model. LMMs are somewhat robust to deviations from normality. 
# Instead, use non-parametric bootstrapping to generate confidence intervals for the fixed effects 
                                            # (does not rely on normality assumption in residuals)

summary(Decomp_model)

##### Non-parametric bootstrapping ----
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

r_squared <- r.squaredGLMM(Decomp_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)







#### SIR ----

hist(final_data_year$CO2CperHourperg_2023)

SIR_model <- lmer(CO2CperHourperg_2023 ~ CO2CperHourperg_2021 + Treatment + Population_Origin * Transplant_Temp + 
                  (1 | Site), 
                  data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(SIR_model)
simulation_output <- simulateResiduals(fittedModel = SIR_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(CO2CperHourperg_2023 ~ CO2CperHourperg_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

# LMM for SIR_model is the final model. The heteroscedasticity is not severe, and there is no indication of overdispersion
# Use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(SIR_model)

##### Non-parametric bootstrapping ----
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

r_squared <- r.squaredGLMM(SIR_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)


















#### SORU Biomass ----

hist(final_data_year$SORU_Biomass_2023)

SORU_biomass_model <- lmer(SORU_Biomass_2023 ~ SORU_Biomass_2021 + Treatment + Population_Origin * Transplant_Temp + 
                           (1 | Site), 
                           data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(SORU_biomass_model)
simulation_output <- simulateResiduals(fittedModel = SORU_biomass_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(SORU_Biomass_2023 ~ SORU_Biomass_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 


# Use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(SORU_biomass_model)

##### Non-parametric bootstrapping ----
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

r_squared <- r.squaredGLMM(SORU_biomass_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)




#### POPRC_Biomass ----
hist(final_data_year$POPRC_Biomass_2023)

POPRC_biomass_model <- lmer(POPRC_Biomass_2023 ~ POPRC_Biomass_2021 + Treatment + Population_Origin * Transplant_Temp + 
                            (1 | Site), 
                            data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(POPRC_biomass_model)
simulation_output <- simulateResiduals(fittedModel = POPRC_biomass_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(POPRC_Biomass_2023 ~ POPRC_Biomass_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 



# Diagnostics look great; still use non-parametric bootstrapping to generate confidence intervals for the fixed effects 

summary(POPRC_biomass_model)

##### Non-parametric bootstrapping ----
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

r_squared <- r.squaredGLMM(POPRC_biomass_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)











#### N-min rate ----
hist(final_data_year$`Overall mineralization rate_2023`)
Nmin_model <- lmer(`Overall mineralization rate_2023` ~ `Overall mineralization rate_2021` + Treatment + Population_Origin * Transplant_Temp + 
                     (1 | Site), 
                   data = final_data_year)
plot(Nmin_model)

# 1. Diagnostics - DHARMa package

simulation_output <- simulateResiduals(fittedModel = Nmin_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) 

# 2. Multicollinearity Assessment
fixed_model <- lm(`Overall mineralization rate_2023` ~ `Overall mineralization rate_2021` + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values)



# Diagnostics look good, but we will still use non-parametric bootstrapping 

summary(Nmin_model)

##### Non-parametric bootstrapping ----
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

r_squared <- r.squaredGLMM(Nmin_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)











#### Soil %C ----
SoilC_model <- lmer(PercentC_SOIL_2023 ~ PercentC_SOIL_2021 + Treatment + Population_Origin * Transplant_Temp + 
                    (1 | Site), 
                    data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(SoilC_model)
simulation_output <- simulateResiduals(fittedModel = SoilC_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PercentC_SOIL_2023 ~ PercentC_SOIL_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values)



summary(SoilC_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(SoilC_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(SoilC_model))
fixed_effects <- fixef(SoilC_model)
SoilC_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(SoilC_effect_summary) # interpret


##### Variance partitioning ----

r_squared <- r.squaredGLMM(SoilC_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)




#### Soil %N ----
SoilN_model <- lmer(PercentN_SOIL_2023 ~ PercentN_SOIL_2021 + Treatment + Population_Origin * Transplant_Temp + 
                    (1 | Site), 
                    data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(SoilN_model)
simulation_output <- simulateResiduals(fittedModel = SoilN_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PercentN_SOIL_2023 ~ PercentN_SOIL_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

summary(SoilN_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(SoilN_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(SoilN_model))
fixed_effects <- fixef(SoilN_model)
SoilN_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(SoilN_effect_summary) # interpret

##### Variance partitioning ----

r_squared <- r.squaredGLMM(SoilN_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)










#### Litter %N ----
LitterN_model <- lmer(PercentN_LITTER_2023 ~ PercentN_LITTER_2021 + Treatment + Population_Origin * Transplant_Temp + 
                      (1 | Site), 
                      data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(LitterN_model)
simulation_output <- simulateResiduals(fittedModel = LitterN_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PercentN_LITTER_2023 ~ PercentN_LITTER_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

summary(LitterN_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(LitterN_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(LitterN_model))
fixed_effects <- fixef(LitterN_model)
LitterN_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(LitterN_effect_summary) # interpret

##### Variance partitioning ----

r_squared <- r.squaredGLMM(LitterN_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)





#### POPRC %N ----

POPRCN_model <- lmer(PercentN_POPRC_2023 ~ PercentN_POPRC_2021 + Treatment + Population_Origin * Transplant_Temp + 
                     (1 | Site), 
                     data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(POPRCN_model)
simulation_output <- simulateResiduals(fittedModel = POPRCN_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PercentN_POPRC_2023 ~ PercentN_POPRC_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

summary(POPRCN_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(POPRCN_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(POPRCN_model))
fixed_effects <- fixef(POPRCN_model)
POPRCN_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(POPRCN_effect_summary) # interpret

##### Variance partitioning ----

r_squared <- r.squaredGLMM(POPRCN_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)




#### SORU %N ----
SORUN_model <- lmer(PercentN_SORU_2023 ~ PercentN_SORU_2021 + Treatment + Population_Origin * Transplant_Temp + 
                    (1 | Site), 
                    data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(SORUN_model)
simulation_output <- simulateResiduals(fittedModel = SORUN_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PercentN_SORU_2023 ~ PercentN_SORU_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 



summary(SORUN_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(SORUN_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(SORUN_model))
fixed_effects <- fixef(SORUN_model)
SORUN_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(SORUN_effect_summary) # interpret

##### Variance partitioning ----

r_squared <- r.squaredGLMM(SORUN_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)










#### Plant Richness ----
PlantRichness_model <- lmer(PlantRichness_2023 ~ PlantRichness_2021 + Treatment + Population_Origin * Transplant_Temp + 
                            (1 | Site), 
                            data = final_data_year)

# 1. Diagnostics - DHARMa package
plot(PlantRichness_model)
simulation_output <- simulateResiduals(fittedModel = PlantRichness_model, n = 1000) 
plot(simulation_output) 
testDispersion(simulation_output) 
testZeroInflation(simulation_output)

# 2. Multicollinearity Assessment
fixed_model <- lm(PlantRichness_2023 ~ PlantRichness_2021 + Treatment + Population_Origin * Transplant_Temp, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

summary(PlantRichness_model)

##### Non-parametric bootstrapping ----
boot_model <- bootMer(PlantRichness_model, FUN = fixef, nsim = 1000)
boot_replicates <- boot_model$t
conf_intervals <- list()
for (i in 1:ncol(boot_replicates)) {
  # Calculate the 2.5th and 97.5th percentiles for the ith fixed effect
  conf_intervals[[i]] <- quantile(boot_replicates[, i], probs = c(0.025, 0.975))
}
conf_intervals_df <- do.call(rbind, conf_intervals)
colnames(conf_intervals_df) <- c("Lower_CI", "Upper_CI")
rownames(conf_intervals_df) <- names(fixef(PlantRichness_model))
fixed_effects <- fixef(PlantRichness_model)
PlantRichness_effect_summary <- data.frame(
  Parameter = names(fixed_effects),
  Estimate = fixed_effects,
  Lower_CI = conf_intervals_df[, "Lower_CI"],
  Upper_CI = conf_intervals_df[, "Upper_CI"]
)

print(PlantRichness_effect_summary) # interpret

##### Variance partitioning ----

r_squared <- r.squaredGLMM(PlantRichness_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)






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

##### Variance partitioning ----

r_squared <- r.squaredGLMM(PlantDiversity_model)
random_effects_variance <- r_squared[2] - r_squared[1]


# Print results
cat("Marginal R² (variance explained by fixed effects):", r_squared[1])
cat("Conditional R² (variance explained by fixed + random effects):", r_squared[2])
cat("Variance explained by random effects:", random_effects_variance)






### SEMs ----

# Load necessary libraries
library(ggplot2)
library(ggpmisc)

# Plot: Change in litter %N and the change in decomposition rate
ggplot(final_data_diff, aes(x = Diff_PercentN_LITTER, y = Diff_MassLoss, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Litter %N vs Change in Decomposition Rate",
       x = "Change in Litter %N",
       y = "Change in Decomposition Rate") +
  facet_wrap(~ Site)

# Plot: Change in litter %N and the change in SIR
ggplot(final_data_diff, aes(x = Diff_PercentN_LITTER, y = Diff_CO2CperHourperg, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Litter %N vs Change in SIR",
       x = "Change in Litter %N",
       y = "Change in SIR") +
  facet_wrap(~ Site)

# Plot: Change in soil %N and the change in SIR
ggplot(final_data_diff, aes(x = Diff_PercentN_SOIL, y = Diff_CO2CperHourperg, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Soil %N vs Change in SIR",
       x = "Change in Soil %N",
       y = "Change in SIR") +
  facet_wrap(~ Site)

# Plot: Change in plant diversity and the change in N-min
ggplot(final_data_diff, aes(x = Diff_PlantDiversity, y = Diff_Overall_mineralization_rate, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Plant Diversity vs Change in N-min",
       x = "Change in Plant Diversity",
       y = "Change in N-min") +
  facet_wrap(~ Site)

# Plot: Change in soil %C and change in decomposition
ggplot(final_data_diff, aes(x = Diff_PercentC_SOIL, y = Diff_MassLoss, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Soil %C vs Change in Decomposition",
       x = "Change in Soil %C",
       y = "Change in Decomposition") +
  facet_wrap(~ Site)

# Plot: Change in soil %N and change in plant diversity
ggplot(final_data_diff, aes(x = Diff_PercentC_SOIL, y = Diff_PlantDiversity, color = Site)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, parse = TRUE, label.x.npc = "right", label.y.npc = "top") +
  labs(title = "Change in Soil %C vs Change in Plant Diversity",
       x = "Change in Soil %C",
       y = "Change in Plant Diversity") +
  facet_wrap(~ Site)
