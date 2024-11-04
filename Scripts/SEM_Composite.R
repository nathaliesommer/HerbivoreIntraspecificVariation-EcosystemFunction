# SEMs - Composite Variable

# Treatment_PopType: Levels (Herbivore-Physiological; Herbivore-Behavioral; Vegetation)

### OVERVIEW: DAG ----

library(ggdag)
library(ggplot2)

# Define paths
# Initial model: notes
dag <- dagify(
  SORUbiomass ~ Treatment_PopType, 
  POPRCbiomass ~ Treatment_PopType,
  SORUN ~ SORUbiomass, #y
  SORUC ~ SORUbiomass, #y
  POPRCC ~ POPRCbiomass, #y
  POPRCN ~ POPRCbiomass, #y
  #LitterC ~ SORUC + POPRCC + MISCC, # removed for singularity + lack of interesting/relevant paths
  LitterN ~ POPRCN + SORUN, #y
  SoilC ~ Treatment_PopType + SORUbiomass + SIR, # removed SOM for collinearity and LitterC for poor fit
  SoilN ~ Treatment_PopType + LitterN, #y
  SIR ~ SoilN + SOM, #y
  #Decomp ~ SIR + LitterN, # removed decomp for poor field proxy
  PlantDiversity ~ SoilC + SOM + SORUbiomass, #y removed SoilN for multicollinearity
  Nmin ~ SoilN + SIR #y (removed SOM for poor fit)
)

# Define paths
# Final model: clean
dag <- dagify(
  SORU_biomass ~ Treatment_PopType, 
  POPRC_biomass ~ Treatment_PopType,
  SORU_N ~ SORU_biomass,
  SORU_C ~ SORU_biomass,
  POPRC_C ~ POPRC_biomass,
  POPRC_N ~ POPRC_biomass,
  Litter_N ~ POPRC_N + SORU_N,
  Soil_C ~ Treatment_PopType + SORU_biomass + SIR, 
  Soil_N ~ Treatment_PopType + LitterN,
  SIR ~ Soil_N + SOM,
  PlantDiversity ~ Soil_C + SOM + SORU_biomass, 
  Nmin ~ Soil_N + SIR 
)

ggdag(dag, text = FALSE) + 
  geom_dag_node(data = . %>% filter(name == "Treatment_PopType"), color = "darkgreen", size = 28) +  
  geom_dag_node(data = . %>% filter(!name == "Treatment_PopType"), size = 28, color = "darkgray") +  
  geom_dag_edges_link(arrow = arrow(length = unit(0.1, "inches")), edge_width = 0.8) +
  geom_dag_text(size = 2.5) + 
  theme_void() +  
  theme(legend.position = "none")

### Scripts to run, objects needed ----

# Decomposition.R: combined_decomp
# SIR.R: SIR_final_data
# N-min.R: N_min_full_data
# Veg.R: functional_groups_wide, combined_diversity_long
# SoilPlantCN.R: CNdata
# LOI.R: SOM_data 

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
  rename(Overall_mineralization_rate = `Overall mineralization rate`) %>% 
  dplyr::select(Sample_ID, Year, Population, Site, Overall_mineralization_rate)

veg_biomass_data <- functional_groups_wide %>%
  rename(Sample_ID = Cage.ID) %>%
  dplyr::select(Sample_ID, Year, Population, Treatment, Transplant, Site, SORU_Biomass, POPRC_Biomass, MISC_Biomass)

diversity_data <- combined_diversity_long %>%
  rename(Sample_ID = Cage.ID) %>%
  mutate(Population = substr(Sample_ID, 1, 2)) %>%  # Add Population column
  dplyr::select(Sample_ID, Year, Population, PlantDiversity)

cn_data <- CNdata

som_data <- SOM_data %>%
  dplyr::select(Sample_ID, average_SOM) %>%
  mutate(Year = 2021)


# Define potential response variables
response_vars <- c(
  "MassLoss",
  "CO2CperHourperg",
  "Overall_mineralization_rate",
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
  "PlantDiversity",
  "average_SOM"
)

# Combine all data frames
combined_data <- decomp_data %>%
  mutate(Year = as.character(Year)) %>%
  left_join(sir_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(n_min_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(veg_biomass_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(diversity_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(cn_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year", "Population")) %>%
  left_join(som_data %>% mutate(Year = as.character(Year)), by = c("Sample_ID", "Year")) %>%
  # Ensure Site column is present
  mutate(Site = coalesce(Site.x, Site.y, Site)) %>%
  # Create a composite Treatment_PopType variable
  mutate(Treatment_PopType = case_when(
    Treatment == "Vegetation" ~ Treatment,
    Population %in% c("YF", "DC", "SC") ~ "Herbivore-Physiological",
    Population %in% c("FN", "UP") ~ "Herbivore-Behavioral"))


# Create final_data_year with variables spread by Year
final_data_year <- combined_data %>%
  pivot_wider(names_from = Year, values_from = all_of(response_vars), names_sep = "_")


final_data_year <- final_data_year %>% 
  select(-Site.y, -Site.x, -average_SOM_2023)


# Relevel for model interpretation
final_data_year$Treatment_PopType <- as.factor(final_data_year$Treatment_PopType)
final_data_year$Treatment_PopType <- relevel(final_data_year$Treatment_PopType, ref = "Vegetation")


# Drop DCS_H3; FNN_H8; SCH_V2 - no data for any metric in 2023
final_data_year <- final_data_year %>% 
  filter(!Sample_ID %in% c("DCS_H3", "FNN_H8", "SCH_V2"))

final_data_year %>% group_by(Sample_ID) %>%
  filter(n() > 1) %>%
  ungroup() %>% print()


# Package list
library(lme4)
library(DHARMa)
library(car)
library(boot)

### Mixed-effects models ----

#### SORU biomass ----

SORU_model <- lmer(SORU_Biomass_2023 ~ 
                     Treatment_PopType +
                     SORU_Biomass_2021 + 
                     (1 | Site), 
                   data = final_data_year)
summary(SORU_model)

##### Check assumptions ----
fixed_model <- lm(SORU_Biomass_2023 ~ 
                    Treatment_PopType + 
                    SORU_Biomass_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SORU_model) # not great
simulation_output <- simulateResiduals(fittedModel = SORU_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # fine

##### Random effects model comparison ---- 
SORU_model_noRE <- lm(SORU_Biomass_2023 ~ 
                        Treatment_PopType +
                        SORU_Biomass_2021,
                      data = final_data_year)

anova(SORU_model, SORU_model_noRE)
VarCorr(SORU_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(SORU_Biomass_2023 ~ 
                Treatment_PopType + 
                SORU_Biomass_2021 + 
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(SORU_model, parm = names(fixef(SORU_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SORU_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SORU_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SORU_model)),
  Original_Estimate = fixef(SORU_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SORU_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) # good

# Print bootstrapped fixed effects comparison
print(boot_summary) # good








#### POPRC biomass ----
POPRC_model <- lmer(POPRC_Biomass_2023 ~ 
                      Treatment_PopType + 
                      POPRC_Biomass_2021 + 
                      (1 | Site), 
                    data = final_data_year)
summary(POPRC_model)

##### Check assumptions ----
fixed_model <- lm(POPRC_Biomass_2023 ~ 
                    Treatment_PopType + 
                    POPRC_Biomass_2021, 
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(POPRC_model)
simulation_output <- simulateResiduals(fittedModel = POPRC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

##### Random effects model comparison ---- 
POPRC_model_noRE <- lm(POPRC_Biomass_2023 ~ 
                         Treatment_PopType +
                         POPRC_Biomass_2021,
                       data = final_data_year)

anova(POPRC_model, POPRC_model_noRE)
VarCorr(POPRC_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(POPRC_Biomass_2023 ~ 
                Treatment_PopType + 
                POPRC_Biomass_2021 +  
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(POPRC_model, parm = names(fixef(POPRC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(POPRC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(POPRC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(POPRC_model)),
  Original_Estimate = fixef(POPRC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(POPRC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) # good

# Print bootstrapped fixed effects comparison
print(boot_summary) # good







#### SORU %N ----
SORUN_model <- lmer(PercentN_SORU_2023 ~ 
                      SORU_Biomass_2023 +
                      PercentN_SORU_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SORUN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_SORU_2023 ~ 
                    SORU_Biomass_2023 + 
                    PercentN_SORU_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SORUN_model) # okay
simulation_output <- simulateResiduals(fittedModel = SORUN_model)
plot(simulation_output) # not good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # zero-inflated

##### Random effects model comparison ---- 
SORUN_model_noRE <- lm(PercentN_SORU_2023 ~ 
                         SORU_Biomass_2023 + 
                         PercentN_SORU_2021,
                       data = final_data_year)

anova(SORUN_model, SORUN_model_noRE)
VarCorr(POPRC_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentN_SORU_2023 ~ 
                SORU_Biomass_2023 + 
                PercentN_SORU_2021 + 
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)



# 1. Extract original confidence intervals
original_ci <- confint(SORUN_model, parm = names(fixef(SORUN_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SORUN_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SORUN_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SORUN_model)),
  Original_Estimate = fixef(SORUN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SORUN_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) # good

# Print bootstrapped fixed effects comparison
print(boot_summary) # good








#### SORU %C ----
SORUC_model <- lmer(PercentC_SORU_2023 ~ 
                      SORU_Biomass_2023 +
                      PercentC_SORU_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SORUC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_SORU_2023 ~ 
                    SORU_Biomass_2023 + 
                    PercentC_SORU_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SORUC_model) # okay
simulation_output <- simulateResiduals(fittedModel = SORUC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # zero-inflated


##### Random effects model comparison ---- 
SORUC_model_noRE <- lm(PercentC_SORU_2023 ~ 
                         SORU_Biomass_2023 + 
                         PercentC_SORU_2021,
                       data = final_data_year)

anova(SORUC_model, SORUC_model_noRE)
VarCorr(SORUC_model)
# keep random effects


##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentC_SORU_2023 ~ 
                SORU_Biomass_2023 + 
                PercentC_SORU_2021 + 
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(SORUC_model, parm = names(fixef(SORUC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SORUC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SORUC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SORUC_model)),
  Original_Estimate = fixef(SORUC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SORUC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)









#### POPRC %C ----
POPRCC_model <- lmer(PercentC_POPRC_2023 ~ 
                       POPRC_Biomass_2023 +
                       PercentC_POPRC_2021 +
                       (1 | Site), 
                     data = final_data_year)
summary(POPRCC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_POPRC_2023 ~ 
                    POPRC_Biomass_2023 + 
                    PercentC_POPRC_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(POPRCC_model) # not great
simulation_output <- simulateResiduals(fittedModel = POPRCC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) #zero-inflated

##### Random effects model comparison ---- 
POPRCC_model_noRE <- lm(PercentC_POPRC_2023 ~ 
                          POPRC_Biomass_2023 + 
                          PercentC_POPRC_2021,
                        data = final_data_year)

anova(POPRCC_model, POPRCC_model_noRE)
VarCorr(POPRCC_model)
# ! drop random effects

POPRCC_model <- lm(PercentC_POPRC_2023 ~ 
                     POPRC_Biomass_2023 +
                     PercentC_POPRC_2021,
                   data = final_data_year)
summary(POPRCC_model)

##### Check assumptions (LM) ----
simulation_output <- simulateResiduals(fittedModel = POPRCC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # zero-inflated



##### Non-parametric bootstrapping for lm ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentC_POPRC_2023 ~ 
              POPRC_Biomass_2023 + 
              PercentC_POPRC_2021, 
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}

# Set seed for reproducibility
set.seed(1231)

# Perform bootstrapping with 1000 resamples
boot_model <- boot(final_data_year, boot_fun, R = 1000) 

# 1. Extract original confidence intervals for lm model using `coef()`
original_ci <- confint(POPRCC_model, parm = names(coef(POPRCC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(POPRCC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(POPRCC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates using `coef()` for lm objects
boot_summary <- data.frame(
  Fixed_Effect = names(coef(POPRCC_model)),
  Original_Estimate = coef(POPRCC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(POPRCC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) 

# Print bootstrapped fixed effects comparison
print(boot_summary) 




#### POPRC %N ----
POPRCN_model <- lmer(PercentN_POPRC_2023 ~ 
                       POPRC_Biomass_2023 +
                       PercentN_POPRC_2021 +
                       (1 | Site), 
                     data = final_data_year)
summary(POPRCN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_POPRC_2023 ~ 
                    POPRC_Biomass_2023 + 
                    PercentN_POPRC_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(POPRCN_model) # good
simulation_output <- simulateResiduals(fittedModel = POPRCN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # zero-inflated


##### Random effects model comparison ---- 
POPRCN_model_noRE <- lm(PercentN_POPRC_2023 ~ 
                          POPRC_Biomass_2023 + 
                          PercentN_POPRC_2021,
                        data = final_data_year)

anova(POPRCN_model, POPRCN_model_noRE)
VarCorr(POPRCN_model)
# ! drop random effects


POPRCN_model <- POPRCN_model_noRE

summary(POPRCN_model)

##### Check assumptions (LM) ----
plot(POPRCN_model) # good
simulation_output <- simulateResiduals(fittedModel = POPRCN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # zero-inflated

##### Non-parametric bootstrapping for lm ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentN_POPRC_2023 ~ 
              POPRC_Biomass_2023 + 
              PercentN_POPRC_2021,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}

# Set seed for reproducibility
set.seed(1231)

# Perform bootstrapping with 1000 resamples
boot_model <- boot(final_data_year, boot_fun, R = 1000) 

# 1. Extract original confidence intervals for lm model using `coef()`
original_ci <- confint(POPRCN_model, parm = names(coef(POPRCN_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(POPRCN_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(POPRCN_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates using `coef()` for lm objects
boot_summary <- data.frame(
  Fixed_Effect = names(coef(POPRCN_model)),
  Original_Estimate = coef(POPRCN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(POPRCN_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) 

# Print bootstrapped fixed effects comparison
print(boot_summary) 











#### SoilN ----

SoilN_model <- lmer(PercentN_SOIL_2023 ~ 
                      Treatment_PopType + 
                      PercentN_LITTER_2023 + 
                      PercentN_SOIL_2021 + 
                      (1 | Site), 
                    data = final_data_year)
summary(SoilN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_SOIL_2023 ~ 
                    Treatment_PopType + 
                    PercentN_LITTER_2023 + 
                    PercentN_SOIL_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SoilN_model) # bad
simulation_output <- simulateResiduals(fittedModel = SoilN_model)
plot(simulation_output) # bad
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good


##### Random effects model comparison ---- 
SoilN_model_noRE <- lm(PercentN_SOIL_2023 ~ 
                         Treatment_PopType + 
                         PercentN_LITTER_2023 + 
                         PercentN_SOIL_2021,
                       data = final_data_year)

anova(SoilN_model, SoilN_model_noRE)
VarCorr(SoilN_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentN_SOIL_2023 ~ 
                Treatment_PopType + 
                PercentN_LITTER_2023 + 
                PercentN_SOIL_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000) 


# 1. Extract original confidence intervals
original_ci <- confint(SoilN_model, parm = names(fixef(SoilN_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SoilN_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SoilN_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SoilN_model)),
  Original_Estimate = fixef(SoilN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SoilN_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)





#### SoilC ----

SoilC_model <- lmer(PercentC_SOIL_2023 ~ 
                      Treatment_PopType + 
                      SORU_Biomass_2023 +
                      CO2CperHourperg_2023 +
                      PercentC_SOIL_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SoilC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_SOIL_2023 ~ 
                    Treatment_PopType + 
                    SORU_Biomass_2023 +
                    CO2CperHourperg_2023 +
                    PercentC_SOIL_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SoilC_model) # not great
simulation_output <- simulateResiduals(fittedModel = SoilC_model)
plot(simulation_output) # bad
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good


##### Random effects model comparison ---- 
SoilC_model_noRE <- lm(PercentC_SOIL_2023 ~ 
                         Treatment_PopType + 
                         SORU_Biomass_2023 +
                         CO2CperHourperg_2023 +
                         PercentC_SOIL_2021,
                       data = final_data_year)

summary(SoilC_model_noRE)

anova(SoilC_model, SoilC_model_noRE)
VarCorr(SoilC_model)
# keep


##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentC_SOIL_2023 ~ 
                Treatment_PopType + 
                SORU_Biomass_2023 +
                CO2CperHourperg_2023 +
                PercentC_SOIL_2021 + 
                (1 | Site),
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000) 


# 1. Extract original confidence intervals
original_ci <- confint(SoilC_model, parm = names(fixef(SoilC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SoilC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SoilC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SoilC_model)),
  Original_Estimate = fixef(SoilC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SoilC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)




#### SIR ----
SIR_model <- lmer(CO2CperHourperg_2023 ~ 
                    PercentN_SOIL_2023 + 
                    average_SOM_2021 +
                    CO2CperHourperg_2021 +
                    (1 | Site), 
                  data = final_data_year)
summary(SIR_model)

##### Check assumptions ----
fixed_model <- lm(CO2CperHourperg_2023 ~ 
                    PercentN_SOIL_2023 + 
                    average_SOM_2021 +
                    CO2CperHourperg_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SIR_model) # okay
simulation_output <- simulateResiduals(fittedModel = SIR_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good 
testZeroInflation(simulation_output) # good


##### Random effects model comparison ---- 
SIR_model_noRE <- lm(CO2CperHourperg_2023 ~ 
                       PercentN_SOIL_2023 + 
                       average_SOM_2021 +
                       CO2CperHourperg_2021,
                     data = final_data_year)

anova(SIR_model, SIR_model_noRE)
VarCorr(SIR_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(CO2CperHourperg_2023 ~ 
                PercentN_SOIL_2023 + 
                average_SOM_2021 +
                CO2CperHourperg_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(SIR_model, parm = names(fixef(SIR_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(SIR_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(SIR_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(SIR_model)),
  Original_Estimate = fixef(SIR_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SIR_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)



#### REMOVE: Decomposition ----
Decomposition_model <- lmer(MassLoss_2023 ~ 
                              CO2CperHourperg_2023 + 
                              PercentN_LITTER_2023 +
                              MassLoss_2021 +
                              (1 | Site), 
                            data = final_data_year)
summary(Decomposition_model)


##### Check assumptions ----
fixed_model <- lm(MassLoss_2023 ~ 
                    CO2CperHourperg_2023 + 
                    PercentN_LITTER_2023 +
                    MassLoss_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(Decomposition_model)
simulation_output <- simulateResiduals(fittedModel = Decomposition_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good


##### Random effects model comparison ---- 
Decomposition_model_noRE <- lm(MassLoss_2023 ~ 
                                 CO2CperHourperg_2023 + 
                                 PercentN_LITTER_2023 +
                                 MassLoss_2021,
                               data = final_data_year)

anova(Decomposition_model, Decomposition_model_noRE)
VarCorr(Decomposition_model)
# ! drop random effects

Decomposition_model <- Decomposition_model_noRE

summary(Decomposition_model)

##### Check assumptions (LM) ----
plot(Decomposition_model)
simulation_output <- simulateResiduals(fittedModel = Decomposition_model)
plot(simulation_output)
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

##### Non-parametric bootstrapping for lm ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(MassLoss_2023 ~ 
              CO2CperHourperg_2023 + 
              PercentN_LITTER_2023 +
              MassLoss_2021,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}

# Set seed for reproducibility
set.seed(1231)

# Perform bootstrapping with 1000 resamples
boot_model <- boot(final_data_year, boot_fun, R = 1000) 

# 1. Extract original confidence intervals for lm model using `coef()`
original_ci <- confint(Decomposition_model, parm = names(coef(Decomposition_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(Decomposition_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(Decomposition_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates using `coef()` for lm objects
boot_summary <- data.frame(
  Fixed_Effect = names(coef(Decomposition_model)),
  Original_Estimate = coef(Decomposition_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(Decomposition_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) 

# Print bootstrapped fixed effects comparison
print(boot_summary) 







#### PlantDiversity ----

PlantDiversity_model <- lmer(PlantDiversity_2023 ~ 
                               PercentC_SOIL_2023 + 
                               SORU_Biomass_2023 + 
                               average_SOM_2021 +
                               PlantDiversity_2021 +
                               (1 | Site), 
                             data = final_data_year)
summary(PlantDiversity_model)


##### Check assumptions ----
fixed_model <- lm(PlantDiversity_2023 ~ 
                    PercentC_SOIL_2023 + 
                    SORU_Biomass_2023 +
                    average_SOM_2021 +
                    PlantDiversity_2021  ,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) 

plot(PlantDiversity_model)
simulation_output <- simulateResiduals(fittedModel = PlantDiversity_model)
plot(simulation_output) # good
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # fine

# given diagnostics, bootstrapping is NOT necessary



##### Random effects model comparison ---- 
PlantDiversity_model_noRE <- lm(PlantDiversity_2023 ~ 
                                  PercentC_SOIL_2023 + 
                                  SORU_Biomass_2023 +
                                  average_SOM_2021 +
                                  PlantDiversity_2021  ,
                                data = final_data_year)

anova(PlantDiversity_model, PlantDiversity_model_noRE)
VarCorr(PlantDiversity_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PlantDiversity_2023 ~ 
                PercentC_SOIL_2023 + 
                average_SOM_2021 +
                SORU_Biomass_2023 + 
                PlantDiversity_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)




# 1. Extract original confidence intervals
original_ci <- confint(PlantDiversity_model, parm = names(fixef(PlantDiversity_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(PlantDiversity_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(PlantDiversity_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(PlantDiversity_model)),
  Original_Estimate = fixef(PlantDiversity_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(PlantDiversity_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)





#### LitterN  ----
LitterN_model <- lmer(PercentN_LITTER_2023 ~ 
                        PercentN_POPRC_2023 + 
                        PercentN_SORU_2023 + 
                        PercentN_LITTER_2021 + 
                        (1 | Site), 
                      data = final_data_year)
summary(LitterN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_LITTER_2023 ~ 
                    PercentN_POPRC_2023 + 
                    PercentN_SORU_2023 + 
                    PercentN_LITTER_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(LitterN_model)
simulation_output <- simulateResiduals(fittedModel = LitterN_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # zero-inflated

##### Random effects model comparison ---- 
LitterN_model_noRE <- lm(PercentN_LITTER_2023 ~ 
                           PercentN_POPRC_2023 + 
                           PercentN_SORU_2023 + 
                           PercentN_LITTER_2021,
                         data = final_data_year)

anova(LitterN_model, LitterN_model_noRE)
VarCorr(LitterN_model)
# keep random effects


##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentN_LITTER_2023 ~ 
                PercentN_POPRC_2023 + 
                PercentN_SORU_2023 + 
                PercentN_LITTER_2021 + 
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(LitterN_model, parm = names(fixef(LitterN_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(LitterN_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(LitterN_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(LitterN_model)),
  Original_Estimate = fixef(LitterN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(LitterN_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)



#### REMOVE: LitterC ----
LitterC_model <- lmer(PercentC_LITTER_2023 ~ 
                        PercentC_MISC_2023 + 
                        PercentC_LITTER_2021 + 
                        (1 | Site), 
                      data = final_data_year)
summary(LitterC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_LITTER_2023 ~ 
                    PercentC_MISC_2023 + 
                    PercentC_LITTER_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(LitterC_model)
simulation_output <- simulateResiduals(fittedModel = LitterC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

# given the diagnostics, non-parametric boostrapping is necessary


##### Random effects model comparison ---- 
LitterC_model_noRE <- lm(PercentC_LITTER_2023 ~ 
                           PercentC_MISC_2023 + 
                           PercentC_LITTER_2021,
                         data = final_data_year)

anova(LitterC_model, LitterC_model_noRE)
VarCorr(LitterC_model)
# keep random effects


##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentC_LITTER_2023 ~ 
                PercentC_MISC_2023 + 
                PercentC_LITTER_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)


# 1. Extract original confidence intervals
original_ci <- confint(LitterC_model, parm = names(fixef(LitterC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(LitterC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(LitterC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(LitterC_model)),
  Original_Estimate = fixef(LitterC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(LitterC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)




#### N-min ----
Nmin_model <- lmer(Overall_mineralization_rate_2023 ~ 
                     PercentN_SOIL_2023 + 
                     CO2CperHourperg_2023 + 
                     Overall_mineralization_rate_2021 +
                     (1 | Site), 
                   data = final_data_year)
summary(Nmin_model)

##### Check assumptions ----
fixed_model <- lm(Overall_mineralization_rate_2023 ~ 
                    PercentN_SOIL_2023 + 
                    CO2CperHourperg_2023 + 
                    Overall_mineralization_rate_2021,
                  data = final_data_year)
vif_values <- vif(fixed_model)
print(vif_values) # good

plot(Nmin_model)
simulation_output <- simulateResiduals(fittedModel = Nmin_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good


##### Random effects model comparison ---- 
Nmin_model_noRE <- lm(Overall_mineralization_rate_2023 ~ 
                        PercentN_SOIL_2023 + 
                        CO2CperHourperg_2023 + 
                        Overall_mineralization_rate_2021,
                      data = final_data_year)

anova(Nmin_model, Nmin_model_noRE)
VarCorr(Nmin_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(Overall_mineralization_rate_2023 ~ 
                PercentN_SOIL_2023 + 
                CO2CperHourperg_2023 + 
                Overall_mineralization_rate_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000)

# 1. Extract original confidence intervals
original_ci <- confint(Nmin_model, parm = names(fixef(Nmin_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(fixef(Nmin_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(fixef(Nmin_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

# 4. Error handling: Ensure that both tables have the same number of rows
if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}

# 5. Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(fixef(Nmin_model)),
  Original_Estimate = fixef(Nmin_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(Nmin_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison)

# Print bootstrapped fixed effects comparison
print(boot_summary)


### SEM-----
library(piecewiseSEM)

# Check individual models
summary(SORU_model)
summary(POPRC_model)
summary(SORUN_model)
summary(SORUC_model)
summary(POPRCC_model)
summary(POPRCN_model)
#summary(LitterC_model)
summary(SoilC_model)
summary(SoilN_model)
summary(SIR_model)
#summary(Decomposition_model)
summary(PlantDiversity_model)
summary(LitterN_model)
summary(Nmin_model)


sem_model_original <- psem(
  SORU_model,
  POPRC_model,
  SORUN_model,
  SORUC_model,
  POPRCC_model,
  POPRCN_model,
  SoilC_model,
  SoilN_model,
  SIR_model,
  PlantDiversity_model,
  LitterN_model,
  Nmin_model,
  data = final_data_year
)

summary(sem_model_original)



