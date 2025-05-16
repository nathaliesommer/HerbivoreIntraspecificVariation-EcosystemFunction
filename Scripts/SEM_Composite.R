# SEMs - Composite Variable

### OVERVIEW: DAG ----

library(ggdag)
library(ggplot2)
library(dplyr)


dag_separate <- dagify(
  Reactor ~ Resistor,
  Resistor ~ Reactor,
  SORU_biomass ~ Reactor + Resistor + Vegetation,
  POPRC_biomass ~ Reactor + Resistor + Vegetation,
  POPRC_biomass ~ SOM,
  SORU_N ~ SORU_biomass,
  SORU_C ~ SORU_biomass,
  POPRC_C ~ POPRC_biomass,
  POPRC_C ~ SORU_biomass,
  POPRC_N ~ POPRC_biomass,
  Litter_N ~ POPRC_N,
  Litter_N ~ SORU_N,
  Soil_C ~ Reactor + Resistor + Vegetation,
  Soil_C ~ SORU_biomass,
  Soil_C ~ SIR,
  Soil_C ~ SOM,
  Soil_N ~ Reactor + Resistor + Vegetation,
  Soil_N ~ Litter_N,
  SIR ~ Soil_N,
  SIR ~ SOM,
  Plant_Diversity ~ Soil_C,
  Plant_Diversity ~ SOM,
  Plant_Diversity ~ SORU_biomass,
  Plant_Diversity ~ POPRC_biomass,
  N_min ~ Soil_N,
  N_min ~ SIR,
  N_min ~ Reactor + Resistor + Vegetation
)

ggdag(dag_separate, text = FALSE) +
  geom_dag_node(data = . %>% filter(name == "Reactor"), color = "#663366", size = 30) +
  geom_dag_node(data = . %>% filter(name == "Resistor"), color = "#D4AF37", size = 30) +
  geom_dag_node(data = . %>% filter(name == "Vegetation"), color = "#9BA48C", size = 30) +
  geom_dag_node(data = . %>% filter(!name %in% c("Reactor", "Resistor", "Vegetation")), size = 28, color = "gray") +
  geom_dag_edges_link(arrow = arrow(length = unit(0.1, "inches")), edge_width = 0.8) +
  geom_dag_text(size = 2.5) +
  theme_void() +
  theme(legend.position = "none")

### Scripts to run, objects needed ----

# SIR.R: SIR_final_data
# N-min.R: N_min_full_data
# Veg.R: functional_groups_wide, combined_diversity_long
# SoilPlantCN.R: CNdata
# LOI.R: SOM_data 

# Load libraries
library(dplyr)
library(tidyr)

### Join data frames ----

# Prepare data from each script

sir_data <- SIR_final_data %>%
  mutate(Population = substr(Sample_ID, 1, 2)) %>%  # Add Population column
  dplyr::select(Sample_ID, Year, Population, CO2CperHourperg)

n_min_data <- N_min_full_data %>%
  rename(Sample_ID = Sample.ID) %>%
  rename(Overall_mineralization_rate = `Overall mineralization rate`) %>% 
  dplyr::select(Sample_ID, Year, Population, Site, Overall_mineralization_rate)

veg_biomass_data <- functional_groups_wide %>%
  rename(Sample_ID = Cage.ID) %>%
  dplyr::select(Sample_ID, Site, Year, Population, Treatment, Transplant, Site, SORU_Biomass, POPRC_Biomass, MISC_Biomass)

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
combined_data <- veg_biomass_data %>%
  mutate(Year = as.character(Year)) %>%
  # Join with other datasets
  left_join(sir_data %>% mutate(Year = as.character(Year)), 
            by = c("Sample_ID", "Year", "Population")) %>%
  left_join(n_min_data %>% mutate(Year = as.character(Year)),
            by = c("Sample_ID", "Year", "Population")) %>%
  left_join(diversity_data %>% mutate(Year = as.character(Year)), 
            by = c("Sample_ID", "Year", "Population")) %>%
  left_join(cn_data %>% mutate(Year = as.character(Year)), 
            by = c("Sample_ID", "Year", "Population")) %>%
  left_join(som_data %>% mutate(Year = as.character(Year)), 
            by = c("Sample_ID", "Year")) %>%
  # Create a composite Treatment_PopType variable
  mutate(Treatment_PopType = case_when(
    Treatment == "Vegetation" ~ Treatment,
    Population %in% c("YF", "DC", "SC") ~ "Herbivore-Physiological",
    Population %in% c("FN", "UP") ~ "Herbivore-Behavioral"))


# Create final_data_year with variables spread by Year
final_data_year <- combined_data %>%
  pivot_wider(names_from = Year, values_from = all_of(response_vars), names_sep = "_")


final_data_year <- final_data_year %>% 
  dplyr::select(-Site.y, -Site.x, -average_SOM_2023)


# Relevel for model interpretation
final_data_year$Treatment_PopType <- as.factor(final_data_year$Treatment_PopType)
final_data_year$Treatment_PopType <- relevel(final_data_year$Treatment_PopType, ref = "Vegetation")


# Drop DCS_H3; FNN_H8; SCH_V2 - no data for any metric in 2023
final_data_year <- final_data_year %>% 
  filter(!Sample_ID %in% c("DCS_H3", "FNN_H8", "SCH_V2"))

# Check
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
                     PercentN_SOIL_2021 +
                     POPRC_Biomass_2021 +
                     PercentN_SORU_2021 +
                     SORU_Biomass_2021 + 
                     (1 | Site), 
                   data = final_data_year)
summary(SORU_model)

##### Check assumptions ----
fixed_model <- lm(SORU_Biomass_2023 ~ 
                    Treatment_PopType +
                    PercentN_SOIL_2021 +
                    POPRC_Biomass_2021 +
                    PercentN_SORU_2021 +
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
                        PercentN_SOIL_2021 +
                        POPRC_Biomass_2021 +
                        PercentN_SORU_2021 +
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
                PercentN_SOIL_2021 +
                POPRC_Biomass_2021 +
                PercentN_SORU_2021 +
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
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(SORU_model) - apply(boot_model$t, 2, median))

# Print bootstrapped CI comparison
print(ci_comparison) # good

# Print bootstrapped fixed effects comparison
print(boot_summary) # good







#### POPRC biomass ----
POPRC_model <- lmer(POPRC_Biomass_2023 ~ 
                      Treatment_PopType + 
                      POPRC_Biomass_2021 + 
                      average_SOM_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(POPRC_model)

##### Check assumptions ----
fixed_model <- lm(POPRC_Biomass_2023 ~ 
                    Treatment_PopType + 
                    POPRC_Biomass_2021 +
                    average_SOM_2021,
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
                         POPRC_Biomass_2021 +
                         average_SOM_2021,
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
                average_SOM_2021 +  
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
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = fixef(POPRC_model) - apply(boot_model$t, 2, median))

# Print bootstrapped CI comparison
print(ci_comparison) # good

# Print bootstrapped fixed effects comparison
print(boot_summary) # good







#### SORU %N ----
SORUN_model <- lmer(PercentN_SORU_2023 ~ 
                      SORU_Biomass_2023 +
                      PlantDiversity_2021 +
                      PercentN_SORU_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SORUN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_SORU_2023 ~ 
                    SORU_Biomass_2023 +
                    PlantDiversity_2021 +
                    PercentN_SORU_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SORUN_model) # okay
simulation_output <- simulateResiduals(fittedModel = SORUN_model)
plot(simulation_output) # not good, warning message expected here
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # zero-inflated

##### Random effects model comparison ---- 
SORUN_model_noRE <- lm(PercentN_SORU_2023 ~ 
                         SORU_Biomass_2023 +
                         PlantDiversity_2021 +
                         PercentN_SORU_2021,
                       data = final_data_year)

anova(SORUN_model, SORUN_model_noRE)
VarCorr(SORUN_model)
# keep random effects

##### Non-parametric bootstrapping ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lmer(PercentN_SORU_2023 ~ 
                SORU_Biomass_2023 +
                PlantDiversity_2021 +
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
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
  if (!is.null(boot_ci_i$percent)) {
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
                      PlantDiversity_2021 +
                      PercentC_SORU_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SORUC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_SORU_2023 ~ 
                    SORU_Biomass_2023 +
                    PlantDiversity_2021 +
                    PercentC_SORU_2021,
                  data = final_data_year)

vif_values <- vif(fixed_model)
print(vif_values) # good

plot(SORUC_model) # okay
simulation_output <- simulateResiduals(fittedModel = SORUC_model)
plot(simulation_output) # not great, warning expected
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # zero-inflated


##### Random effects model comparison ---- 
SORUC_model_noRE <- lm(PercentC_SORU_2023 ~ 
                         SORU_Biomass_2023 +
                         PlantDiversity_2021 +
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
                PlantDiversity_2021 +
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
                       SORU_Biomass_2023 +
                       PlantDiversity_2021 +
                       Overall_mineralization_rate_2021 +
                       POPRC_Biomass_2021 +
                       PercentC_SOIL_2021 +
                       PercentC_POPRC_2021 +
                       (1 | Site), 
                     data = final_data_year)
summary(POPRCC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_POPRC_2023 ~ 
                    POPRC_Biomass_2023 +
                    SORU_Biomass_2023 +
                    PlantDiversity_2021 +
                    Overall_mineralization_rate_2021 +
                    POPRC_Biomass_2021 +
                    PercentC_SOIL_2021 +
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
                          SORU_Biomass_2023 +
                          PlantDiversity_2021 +
                          Overall_mineralization_rate_2021 +
                          POPRC_Biomass_2021 +
                          PercentC_SOIL_2021 +
                          PercentC_POPRC_2021,
                        data = final_data_year)

anova(POPRCC_model, POPRCC_model_noRE)
VarCorr(POPRCC_model)
# ! drop random effects

POPRCC_model <- POPRCC_model_noRE
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
              SORU_Biomass_2023 +
              PlantDiversity_2021 +
              Overall_mineralization_rate_2021 +
              POPRC_Biomass_2021 +
              PercentC_SOIL_2021 +
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
                       Overall_mineralization_rate_2021 +
                       PercentN_POPRC_2021 +
                       (1 | Site), 
                     data = final_data_year)
# singularity warning; random effects dropped below

summary(POPRCN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_POPRC_2023 ~ 
                    POPRC_Biomass_2023 + 
                    Overall_mineralization_rate_2021 +
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
                          Overall_mineralization_rate_2021 +
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
              Overall_mineralization_rate_2021 +
              PercentN_POPRC_2021,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}

# Set seed for reproducibility
set.seed(1231)

# Perform bootstrapping with 1000 resamples
boot_model <- boot(final_data_year, boot_fun, R = 1000) 

# 1. Extract original confidence intervals for lm model using `coef()`
original_ci <- confint(POPRCN_model, level = 0.95)
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
  tryCatch({
    boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
    if (!is.null(boot_ci_i$percent)) {
      boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4]
      boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
    }
  }, error = function(e) {
    warning(paste("Could not compute CI for coefficient", i))
  })
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
                      average_SOM_2021 +
                      CO2CperHourperg_2021 +
                      (1 | Site), 
                    data = final_data_year)
summary(SoilN_model)

##### Check assumptions ----
fixed_model <- lm(PercentN_SOIL_2023 ~ 
                    Treatment_PopType + 
                    PercentN_LITTER_2023 + 
                    PercentN_SOIL_2021 + 
                    average_SOM_2021 +
                    CO2CperHourperg_2021,
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
                         PercentN_SOIL_2021 + 
                         average_SOM_2021 +
                         CO2CperHourperg_2021, 
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
                average_SOM_2021 +
                CO2CperHourperg_2021 +
                (1 | Site), 
              data = d)
  return(fixef(mod))
}

set.seed(1231)

boot_model <- boot(final_data_year, boot_fun, R = 1000) # warning expected


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
                      CO2CperHourperg_2021 +
                      average_SOM_2021 +
                      PercentC_POPRC_2021 +
                      PercentC_SOIL_2021 + 
                      (1 | Site), 
                    data = final_data_year)
summary(SoilC_model)

##### Check assumptions ----
fixed_model <- lm(PercentC_SOIL_2023 ~ 
                    Treatment_PopType + 
                    SORU_Biomass_2023 +
                    CO2CperHourperg_2023 +
                    CO2CperHourperg_2021 +
                    average_SOM_2021 +
                    PercentC_POPRC_2021 +
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
                         CO2CperHourperg_2021 +
                         PercentC_POPRC_2021 +
                         average_SOM_2021 +
                         PercentC_SOIL_2021,
                       data = final_data_year)

summary(SoilC_model_noRE)

anova(SoilC_model, SoilC_model_noRE)
VarCorr(SoilC_model)
# drop random effects

SoilC_model <- SoilC_model_noRE

##### Check assumptions (LM) ----
simulation_output <- simulateResiduals(fittedModel = SoilC_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # zero-inflated



##### Non-parametric bootstrapping for lm ----
# Function for bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentC_SOIL_2023 ~ 
              Treatment_PopType + 
              SORU_Biomass_2023 +
              CO2CperHourperg_2023 +
              CO2CperHourperg_2021 +
              PercentC_POPRC_2021 +
              average_SOM_2021 +
              PercentC_SOIL_2021,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}

# Set seed for reproducibility
set.seed(1231)

# Perform bootstrapping with 1000 resamples
boot_model <- boot(final_data_year, boot_fun, R = 1000) 

# 1. Extract original confidence intervals for lm model using `coef()`
original_ci <- confint(SoilC_model, parm = names(coef(SoilC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# 2. Initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(SoilC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# 3. Iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(SoilC_model))) {
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
  Fixed_Effect = names(coef(SoilC_model)),
  Original_Estimate = coef(SoilC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(SoilC_model) - apply(boot_model$t, 2, median)
)

# Print bootstrapped CI comparison
print(ci_comparison) 

# Print bootstrapped fixed effects comparison
print(boot_summary) 



#### SIR ----
SIR_model <- lmer(CO2CperHourperg_2023 ~ 
                    PercentN_SOIL_2023 + 
                    average_SOM_2021 + 
                    PlantDiversity_2021 +
                    PercentN_SOIL_2021 +
                    PercentN_POPRC_2021 +
                    POPRC_Biomass_2021 +
                    CO2CperHourperg_2021 + 
                    PercentN_POPRC_2021 +
                    (1 | Site), 
                  data = final_data_year)
summary(SIR_model)

##### Check assumptions ----
fixed_model <- lm(CO2CperHourperg_2023 ~ 
                    PercentN_SOIL_2023 + 
                    average_SOM_2021 + 
                    PlantDiversity_2021 +
                    PercentN_SOIL_2021 +
                    PercentN_POPRC_2021 +
                    POPRC_Biomass_2021 +
                    CO2CperHourperg_2021 + 
                    PercentN_POPRC_2021,
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
                     PlantDiversity_2021 +
                     PercentN_SOIL_2021 +
                     PercentN_POPRC_2021 +
                     POPRC_Biomass_2021 +
                     CO2CperHourperg_2021 + 
                     PercentN_POPRC_2021,
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
                PlantDiversity_2021 +
                PercentN_SOIL_2021 +
                PercentN_POPRC_2021 +
                POPRC_Biomass_2021 +
                CO2CperHourperg_2021 + 
                PercentN_POPRC_2021 +
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





#### PlantDiversity ----

PlantDiversity_model <- lmer(PlantDiversity_2023 ~ 
                               PercentC_SOIL_2023 + 
                               SORU_Biomass_2023 + 
                               POPRC_Biomass_2023 +
                               PercentN_SORU_2021 +
                               average_SOM_2021 +
                               PlantDiversity_2021 +
                               Overall_mineralization_rate_2021 +
                               PercentN_LITTER_2021 +
                               (1 | Site), 
                             data = final_data_year)
summary(PlantDiversity_model)


##### Check assumptions ----
fixed_model <- lm(PlantDiversity_2023 ~ 
                    PercentC_SOIL_2023 + 
                    SORU_Biomass_2023 +
                    POPRC_Biomass_2023 +
                    PercentN_SORU_2021 +
                    average_SOM_2021 +
                    PlantDiversity_2021 +
                    PercentN_LITTER_2021 +
                    Overall_mineralization_rate_2021,
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
                                  POPRC_Biomass_2023 +
                                  PercentN_SORU_2021 +
                                  average_SOM_2021 +
                                  Overall_mineralization_rate_2021 +
                                  PercentN_LITTER_2021 +
                                  PlantDiversity_2021,
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
                SORU_Biomass_2023 +
                POPRC_Biomass_2023 +
                PercentN_SORU_2021 +
                average_SOM_2021 +
                Overall_mineralization_rate_2021 +
                PercentN_LITTER_2021 +
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


#### N-min ----
Nmin_model <- lmer(Overall_mineralization_rate_2023 ~ 
                     Treatment_PopType +
                     PercentN_SOIL_2023 + 
                     CO2CperHourperg_2023 + 
                     PercentN_POPRC_2021 +
                     Overall_mineralization_rate_2021 +
                     (1 | Site), 
                   data = final_data_year)
summary(Nmin_model)

##### Check assumptions ----
fixed_model <- lm(Overall_mineralization_rate_2023 ~ 
                    Treatment_PopType +
                    PercentN_SOIL_2023 + 
                    CO2CperHourperg_2023 + 
                    PercentN_POPRC_2021 +
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
                        Treatment_PopType +
                        PercentN_SOIL_2023 + 
                        CO2CperHourperg_2023 + 
                        PercentN_POPRC_2021 +
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
                Treatment_PopType +
                PercentN_SOIL_2023 + 
                CO2CperHourperg_2023 + 
                PercentN_POPRC_2021 +
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
library(lme4)

# Check individual models
summary(SORU_model)
summary(POPRC_model)
summary(SORUN_model)
summary(SORUC_model)
summary(POPRCC_model)
summary(POPRCN_model)
summary(SoilC_model)
summary(SoilN_model)
summary(SIR_model)
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

summary_sem_model <- summary(sem_model_original)

options(max.print = 10000)
sink("SEM_Summary_Final.txt") # saved as SEM_Summary_Initial before incorporating tests of directed separation 
print(summary_sem_model)
sink()


### Extract R2 for random effects interpretation----
# Load necessary library
if (!require(MuMIn)) install.packages("MuMIn")
library(MuMIn)

# List of models
model_list <- list(
  SORU_model = SORU_model,
  POPRC_model = POPRC_model,
  SORUN_model = SORUN_model,
  SORUC_model = SORUC_model,
  POPRCC_model = POPRCC_model,
  POPRCN_model = POPRCN_model,
  SoilC_model = SoilC_model,
  SoilN_model = SoilN_model,
  SIR_model = SIR_model,
  PlantDiversity_model = PlantDiversity_model,
  LitterN_model = LitterN_model,
  Nmin_model = Nmin_model
)

# Create a data frame to store the R-squared values
r2_summary <- data.frame(
  Model = character(),
  Marginal_R2 = numeric(),
  Conditional_R2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each model and extract R-squared values
for (model_name in names(model_list)) {
  model <- model_list[[model_name]]
  
     # Extract R-squared values using MuMIn
     r2_values <- r.squaredGLMM(model)
     
     # Add the R-squared values to the summary data frame
     r2_summary <- rbind(r2_summary, data.frame(
       Model = model_name,
       Marginal_R2 = r2_values[1],
       Conditional_R2 = r2_values[2]
    ))
}

# Print the final R-squared summary for debugging
print(r2_summary)

# Define file path
output_file <- "R_squared_summary.txt"

# Write the summary to a text file
write.table(r2_summary, file = output_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)




### Extract standardized coefficients for indirect path assessment ----
coefficients_table <- summary_sem_model$coefficients

coefficients_table$Estimate <- as.numeric(coefficients_table$Estimate) #NAs warning expected

filtered_coefficients <- coefficients_table[coefficients_table$P.Value < 0.05, ] # filter for significance

result <- filtered_coefficients[, c("Response", "Predictor", "Estimate")]



# Step 1: Identify Intermediate Variables
treatment_categories <- c("Treatment_PopType = Herbivore-Behavioral", 
                          "Treatment_PopType = Herbivore-Physiological", 
                          "Treatment_PopType = Vegetation")

# Find direct paths from treatment categories
direct_paths <- result[result$Predictor %in% treatment_categories, ]

# Function to recursively find indirect paths
find_indirect_paths <- function(current_path, current_estimate, current_intermediate) {
  # Find downstream paths from the current intermediate
  downstream <- result[result$Predictor == current_intermediate, ]
  
  # If no downstream paths, return the current path
  if (nrow(downstream) == 0) {
    return(list(current_path))
  }
  
  # List to store all paths
  all_paths <- list()
  
  # Iterate over each downstream path
  for (i in seq_len(nrow(downstream))) {
    next_intermediate <- downstream$Response[i]
    next_estimate <- downstream$Estimate[i]
    
    # Calculate new indirect effect size
    new_estimate <- current_estimate * next_estimate
    
    # Create new path
    new_path <- append(current_path, list(
      list(
        Intermediate = next_intermediate,
        Estimate = next_estimate,
        Indirect_Effect_Size = new_estimate
      )
    ))
    
    # Recursively find further paths
    all_paths <- c(all_paths, find_indirect_paths(new_path, new_estimate, next_intermediate))
  }
  
  return(all_paths)
}

# Step 2: Map Paths for Indirect Influence
# Create a list to store all indirect paths
all_indirect_paths <- list()

# Iterate over each direct path
for (i in seq_len(nrow(direct_paths))) {
  treatment <- direct_paths$Predictor[i]
  intermediate <- direct_paths$Response[i]
  estimate1 <- direct_paths$Estimate[i]
  
  # Initialize the path
  initial_path <- list(
    list(
      Treatment = treatment,
      Intermediate = intermediate,
      Estimate = estimate1,
      Indirect_Effect_Size = estimate1
    )
  )
  
  # Find all indirect paths starting from this direct path
  paths <- find_indirect_paths(initial_path, estimate1, intermediate)
  
  # Add to the list of all indirect paths
  all_indirect_paths <- c(all_indirect_paths, paths)
}

# Step 3: Convert the list of indirect paths to a data frame
indirect_effects_df <- do.call(rbind, lapply(all_indirect_paths, function(path) {
  # Flatten the path into a single row
  data.frame(
    Treatment = path[[1]]$Treatment,
    Intermediates = paste(sapply(path, function(x) x$Intermediate), collapse = " -> "),
    Indirect_Effect_Size = path[[length(path)]]$Indirect_Effect_Size
  )
}))

# Print the indirect effects
print(indirect_effects_df)


# Define file path
output_file <- "Indirect_paths.txt"

# Write the summary to a text file
write.table(indirect_effects_df, file = output_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# Figures ----

library(DiagrammeR)

# Define the mapping for labels
label_mapping <- c(
  "SORU_Biomass_2023" = "Goldenrod biomass",
  "POPRC_Biomass_2023" = "Grass biomass",
  "PercentN_SORU_2023" = "Goldenrod %N",
  "PercentC_SORU_2023" = "Goldenrod %C",
  "PercentN_POPRC_2023" = "Grass %N",
  "PercentC_POPRC_2023" = "Grass %C",
  "PercentC_SOIL_2023" = "Soil %C",
  "PercentN_SOIL_2023" = "Soil %N",
  "CO2CperHourperg_2023" = "SIR",
  "PlantDiversity_2023" = "Plant diversity",
  "PercentN_LITTER_2023" = "Litter %N",
  "Overall_mineralization_rate_2023" = "Nitrogen mineralization",
  "Herbivore" = "Herbivore",
  "Treatment_PopType = Vegetation" = "Vegetation",
  "average_SOM_2021" = "Baseline SOM",
  "SORU_Biomass_2021" = "Baseline goldenrod biomass",
  "POPRC_Biomass_2021" = "Baseline grass biomass",
  "PercentN_SORU_2021" = "Baseline goldenrod %N",
  "PercentC_SORU_2021" = "Baseline goldenrod %C",
  "PercentN_POPRC_2021" = "Baseline grass %N",
  "PercentC_POPRC_2021" = "Baseline grass %C",
  "PercentC_SOIL_2021" = "Baseline soil %C",
  "PercentN_SOIL_2021" = "Baseline soil %N",
  "CO2CperHourperg_2021" = "Baseline SIR",
  "PlantDiversity_2021" = "Baseline plant diversity",
  "PercentN_LITTER_2021" = "Baseline litter %N",
  "Overall_mineralization_rate_2021" = "Baseline nitrogen mineralization"
)

#### Basic Significant Paths ----

# Define the mapping for labels (Plasticity version)
plasticity_label_mapping <- c(
  "SORU_Biomass_2023" = "Goldenrod biomass",
  "POPRC_Biomass_2023" = "Grass biomass",
  "PercentN_SORU_2023" = "Goldenrod %N",
  "PercentC_SORU_2023" = "Goldenrod %C",
  "PercentN_POPRC_2023" = "Grass %N",
  "PercentC_POPRC_2023" = "Grass %C",
  "PercentC_SOIL_2023" = "Soil %C",
  "PercentN_SOIL_2023" = "Soil %N",
  "CO2CperHourperg_2023" = "SIR",
  "PlantDiversity_2023" = "Plant diversity",
  "PercentN_LITTER_2023" = "Litter %N",
  "Overall_mineralization_rate_2023" = "Nitrogen mineralization",
  "Treatment_PopType = Herbivore-Behavioral" = "Reactors",
  "Treatment_PopType = Herbivore-Physiological" = "Resistors",
  "Treatment_PopType = Vegetation" = "Vegetation",
  "average_SOM_2021" = "Baseline SOM",
  "SORU_Biomass_2021" = "Baseline goldenrod biomass",
  "POPRC_Biomass_2021" = "Baseline grass biomass",
  "PercentN_SORU_2021" = "Baseline goldenrod %N",
  "PercentC_SORU_2021" = "Baseline goldenrod %C",
  "PercentN_POPRC_2021" = "Baseline grass %N",
  "PercentC_POPRC_2021" = "Baseline grass %C",
  "PercentC_SOIL_2021" = "Baseline soil %C",
  "PercentN_SOIL_2021" = "Baseline soil %N",
  "CO2CperHourperg_2021" = "Baseline SIR",
  "PlantDiversity_2021" = "Baseline plant diversity",
  "PercentN_LITTER_2021" = "Baseline litter %N",
  "Overall_mineralization_rate_2021" = "Baseline nitrogen mineralization"
)

# Clean up the data frame and apply label mapping
clean_coefficients <- data.frame(
  Predictor = filtered_coefficients$Predictor,
  Response = filtered_coefficients$Response,
  Estimate = filtered_coefficients$Estimate,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Predictor) & !is.na(Response) & !is.na(Estimate)) %>%
  mutate(
    Predictor_Label = plasticity_label_mapping[Predictor],
    Response_Label = plasticity_label_mapping[Response],
    # Calculate penwidth based on absolute effect size
    Penwidth = (abs(Estimate) / max(abs(Estimate))) * 6 + 1  # Scale from 1 to 7
  )

# Create basic dot script
basic_dot_script <- paste0(
  "digraph SEM_PathDiagram { \n",
  "  rankdir=LR;\n",
  "  node [shape=rectangle];\n",
  "  splines=true;\n",
  
  # Add baseline nodes with gray fill
  paste0("  \"", unique(grep("Baseline", clean_coefficients$Predictor_Label, value=TRUE)), 
         "\" [style=filled, fillcolor=gray85];\n", collapse=""),
  "\n",
  
  # Add edges with dotted lines for negative effects and scaled width
  paste0("  \"", clean_coefficients$Predictor_Label, "\" -> \"", 
         clean_coefficients$Response_Label, "\" [label=\" ", 
         round(clean_coefficients$Estimate, 2),
         "\"", 
         ifelse(clean_coefficients$Estimate < 0, ", style=dashed", ""),
         ", penwidth=", round(clean_coefficients$Penwidth, 2),
         "];\n",
         collapse = ""),
  
  "}"
)

# Render the plasticity plot
grViz(basic_dot_script)


#### Herbivory ----

filtered_coefficients_herbivory <- filtered_coefficients %>%
  dplyr::select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  # First group by Response to handle the herbivore types averaging
  group_by(Response) %>%
  mutate(
    Estimate = ifelse(
      Predictor %in% c("Treatment_PopType = Herbivore-Behavioral", "Treatment_PopType = Herbivore-Physiological"),
      mean(Estimate[Predictor %in% c("Treatment_PopType = Herbivore-Behavioral", "Treatment_PopType = Herbivore-Physiological")]),
      Estimate
    ),
    Predictor = ifelse(
      Predictor %in% c("Treatment_PopType = Herbivore-Behavioral", "Treatment_PopType = Herbivore-Physiological"),
      "Herbivore",
      Predictor
    )
  ) %>%
  ungroup() %>%
  # Now filter baseline variables
  group_by(Response) %>%
  filter(
    # Keep non-baseline variables
    !grepl("2021$", Predictor) |
    # Or keep baseline variables only if their absolute effect size is greater than non-baseline effects
    (grepl("2021$", Predictor) & abs(Estimate) >= {
      non_baseline_estimates <- abs(Estimate[!grepl("2021$", Predictor)])
      if(length(non_baseline_estimates) > 0) {
        max(non_baseline_estimates) + 1
      } else {
        -Inf  # If no non-baseline predictors, keep all baseline predictors
      }
    })
  ) %>%
  distinct(Response, Predictor, .keep_all = TRUE)




# Apply label mapping to predictors and responses
filtered_coefficients_herbivory_mapped <- filtered_coefficients_herbivory %>%
  mutate(
    Predictor_Label = label_mapping[Predictor],
    Response_Label = label_mapping[Response]
  )

# Check for any unmapped labels
unmapped_predictors <- filtered_coefficients_herbivory$Predictor[is.na(filtered_coefficients_herbivory_mapped$Predictor_Label)]
unmapped_responses <- filtered_coefficients_herbivory$Response[is.na(filtered_coefficients_herbivory_mapped$Response_Label)]

if(length(unmapped_predictors) > 0){
  warning("Unmapped Predictors found: ", paste(unmapped_predictors, collapse = ", "))
}

if(length(unmapped_responses) > 0){
  warning("Unmapped Responses found: ", paste(unmapped_responses, collapse = ", "))
}

# Optionally, remove rows with unmapped labels
filtered_coefficients_herbivory_mapped <- filtered_coefficients_herbivory_mapped %>%
  filter(!is.na(Predictor_Label) & !is.na(Response_Label))

# Create Display Labels with Correct Spacing for Baseline Nodes
filtered_coefficients_herbivory_mapped <- filtered_coefficients_herbivory_mapped %>%
  mutate(
    Predictor_Display_Label = ifelse(
      grepl("Baseline", Predictor_Label, ignore.case = TRUE),
      sub("(?i)(Baseline)(\\S)", "Baseline \\2", Predictor_Label, perl = TRUE),
      Predictor_Label
    ),
    Response_Display_Label = ifelse(
      grepl("Baseline", Response_Label, ignore.case = TRUE),
      sub("(?i)(Baseline)(\\S)", "Baseline \\2", Response_Label, perl = TRUE),
      Response_Label
    )
  )

# Calculate absolute estimates for scaling
filtered_coefficients_herbivory_mapped <- filtered_coefficients_herbivory_mapped %>%
  mutate(
    Absolute_Estimate = abs(Estimate)
  )

# Define minimum and maximum penwidth
min_penwidth <- 2
max_penwidth <- 7.5

# Calculate scaling factor
max_abs_estimate <- max(filtered_coefficients_herbivory_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_herbivory_mapped <- filtered_coefficients_herbivory_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Function to determine edge color based on Estimate
get_edge_color <- function(estimate) {
  if (estimate < 0) {
    return("red")
  } else {
    return("black")
  }
}

# Function to find all downstream nodes and paths
find_downstream_paths <- function(start_nodes, all_paths) {
  # Initialize a vector to keep track of visited nodes
  visited <- c()
  # Initialize a queue with the start nodes
  queue <- start_nodes
  
  # While there are nodes to process
  while (length(queue) > 0) {
    # Pop the first node from the queue
    current_node <- queue[1]
    queue <- queue[-1]
    
    # Mark the current node as visited
    visited <- unique(c(visited, current_node))
    
    # Find all paths starting from the current node
    downstream_paths <- all_paths %>%
      filter(Predictor_Display_Label == current_node)
    
    # Add the response nodes to the queue if they haven't been visited
    new_nodes <- setdiff(downstream_paths$Response_Display_Label, visited)
    queue <- unique(c(queue, new_nodes))
  }
  
  # Return all paths that involve the visited nodes
  all_paths %>%
    filter(Predictor_Display_Label %in% visited | Response_Display_Label %in% visited)
}

# Filter to include only paths starting from "Herbivore" or "Vegetation" and their downstream paths
filtered_paths <- find_downstream_paths(
  start_nodes = c("Herbivore", "Vegetation"),
  all_paths = filtered_coefficients_herbivory_mapped
)

# Extract unique display labels for nodes involved in these paths
unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

# Initialize the DOT script
dot_script <- "digraph SEM_PathDiagram { 
  rankdir=LR;
  node [shape=rectangle, style=filled, fillcolor=lightblue];
  
"

# Add nodes with shape customization and background color
for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    # Use rectangle shape with dark gray background for baseline nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    # Use default rectangle shape with white background for other nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

# Add edges with labels, colors, and penwidth
for(i in 1:nrow(filtered_paths)){
  predictor <- filtered_paths$Predictor_Display_Label[i]
  response <- filtered_paths$Response_Display_Label[i]
  estimate <- filtered_paths$Estimate[i]
  color <- get_edge_color(estimate)
  penwidth <- filtered_paths$Penwidth[i]
  label <- sprintf("%.2f", estimate)
  
  dot_script <- paste0(dot_script, sprintf("  \"%s\" -> \"%s\" [color=%s, label=\"%s\", penwidth=%.2f];\n",
                                           predictor, response, color, label, penwidth))
}

# Close the DOT script
dot_script <- paste0(dot_script, "}")

# Render the graph using grViz
grViz(dot_script)

#### Plasticity ----

filtered_coefficients_plasticity <- filtered_coefficients %>%
  dplyr::select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  group_by(Response) %>%
  # Subtract the estimate for vegetation from herbivore types
  mutate(
    Vegetation_Estimate = ifelse(
      any(Predictor == "Treatment_PopType = Vegetation"),
      Estimate[Predictor == "Treatment_PopType = Vegetation"],
      0  # Default value if Vegetation is not present
    ),
    Estimate = ifelse(
      Predictor %in% c("Treatment_PopType = Herbivore-Behavioral", "Treatment_PopType = Herbivore-Physiological"),
      Estimate - Vegetation_Estimate,
      Estimate
    ),
    Predictor = case_when(
      Predictor == "Treatment_PopType = Herbivore-Behavioral" ~ "Behaviorally plastic",
      Predictor == "Treatment_PopType = Herbivore-Physiological" ~ "Physiologically plastic",
      TRUE ~ Predictor
    )
  ) %>%
  # Exclude "Treatment_PopType = Vegetation"
  filter(Predictor != "Treatment_PopType = Vegetation") %>%
  # Filter to include "2021" predictors only if their absolute effect size is greater by at least 1
  filter(
    !grepl("2021$", Predictor) | 
      (grepl("2021$", Predictor) & abs(Estimate) >= {
        non_baseline_estimates <- abs(Estimate[!grepl("2021$", Predictor)])
        if (length(non_baseline_estimates) > 0) {
          max(non_baseline_estimates) + 1
        } else {
          Inf  # If no non-baseline predictors, set a high threshold to exclude
        }
      })
  ) %>%
  # Remove duplicate rows after averaging
  distinct()

# Re-define the mapping for labels
label_mapping <- c(
  "SORU_Biomass_2023" = "Goldenrod biomass",
  "POPRC_Biomass_2023" = "Grass biomass",
  "PercentN_SORU_2023" = "Goldenrod %N",
  "PercentC_SORU_2023" = "Goldenrod %C",
  "PercentN_POPRC_2023" = "Grass %N",
  "PercentC_POPRC_2023" = "Grass %C",
  "PercentC_SOIL_2023" = "Soil %C",
  "PercentN_SOIL_2023" = "Soil %N",
  "CO2CperHourperg_2023" = "SIR",
  "PlantDiversity_2023" = "Plant diversity",
  "PercentN_LITTER_2023" = "Litter %N",
  "Overall_mineralization_rate_2023" = "Nitrogen mineralization",
  "Behaviorally plastic" = "Behaviorally plastic",
  "Physiologically plastic" = "Physiologically plastic",
  "Treatment_PopType = Vegetation" = "Vegetation",
  "average_SOM_2021" = "Baseline SOM",
  "SORU_Biomass_2021" = "Baseline goldenrod biomass",
  "POPRC_Biomass_2021" = "Baseline grass biomass",
  "PercentN_SORU_2021" = "Baseline goldenrod %N",
  "PercentC_SORU_2021" = "Baseline goldenrod %C",
  "PercentN_POPRC_2021" = "Baseline grass %N",
  "PercentC_POPRC_2021" = "Baseline grass %C",
  "PercentC_SOIL_2021" = "Baseline soil %C",
  "PercentN_SOIL_2021" = "Baseline soil %N",
  "CO2CperHourperg_2021" = "Baseline SIR",
  "PlantDiversity_2021" = "Baseline plant diversity",
  "PercentN_LITTER_2021" = "Baseline litter %N",
  "Overall_mineralization_rate_2021" = "Baseline nitrogen mineralization"
)

# Apply label mapping to predictors and responses
filtered_coefficients_plasticity_mapped <- filtered_coefficients_plasticity %>%
  mutate(
    Predictor_Label = label_mapping[Predictor],
    Response_Label = label_mapping[Response]
  )

# Check for any unmapped labels
unmapped_predictors <- filtered_coefficients_plasticity$Predictor[is.na(filtered_coefficients_plasticity_mapped$Predictor_Label)]
unmapped_responses <- filtered_coefficients_plasticity$Response[is.na(filtered_coefficients_plasticity_mapped$Response_Label)]

if(length(unmapped_predictors) > 0){
  warning("Unmapped Predictors found: ", paste(unmapped_predictors, collapse = ", "))
}

if(length(unmapped_responses) > 0){
  warning("Unmapped Responses found: ", paste(unmapped_responses, collapse = ", "))
}

# Optionally, remove rows with unmapped labels
filtered_coefficients_plasticity_mapped <- filtered_coefficients_plasticity_mapped %>%
  filter(!is.na(Predictor_Label) & !is.na(Response_Label))

# Create Display Labels with Correct Spacing for Baseline Nodes
filtered_coefficients_plasticity_mapped <- filtered_coefficients_plasticity_mapped %>%
  mutate(
    Predictor_Display_Label = ifelse(
      grepl("Baseline", Predictor_Label, ignore.case = TRUE),
      sub("(?i)(Baseline)(\\S)", "Baseline \\2", Predictor_Label, perl = TRUE),
      Predictor_Label
    ),
    Response_Display_Label = ifelse(
      grepl("Baseline", Response_Label, ignore.case = TRUE),
      sub("(?i)(Baseline)(\\S)", "Baseline \\2", Response_Label, perl = TRUE),
      Response_Label
    )
  )

# Calculate absolute estimates for scaling
filtered_coefficients_plasticity_mapped <- filtered_coefficients_plasticity_mapped %>%
  mutate(
    Absolute_Estimate = abs(Estimate)
  )

# Define minimum and maximum penwidth
min_penwidth <- 2
max_penwidth <- 7.5

# Calculate scaling factor
max_abs_estimate <- max(filtered_coefficients_plasticity_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_plasticity_mapped <- filtered_coefficients_plasticity_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Function to find all downstream nodes and paths
find_downstream_paths <- function(start_nodes, all_paths) {
  # Initialize a vector to keep track of visited nodes
  visited <- c()
  # Initialize a queue with the start nodes
  queue <- start_nodes
  
  # While there are nodes to process
  while (length(queue) > 0) {
    # Pop the first node from the queue
    current_node <- queue[1]
    queue <- queue[-1]
    
    # Mark the current node as visited
    visited <- unique(c(visited, current_node))
    
    # Find all paths starting from the current node
    downstream_paths <- all_paths %>%
      filter(Predictor_Display_Label == current_node)
    
    # Add the response nodes to the queue if they haven't been visited
    new_nodes <- setdiff(downstream_paths$Response_Display_Label, visited)
    queue <- unique(c(queue, new_nodes))
  }
  
  # Return all paths that involve the visited nodes
  all_paths %>%
    filter(Predictor_Display_Label %in% visited | Response_Display_Label %in% visited)
}

# Filter to include only paths starting from "Behaviorally plastic" or "Physiologically plastic" and their downstream paths
filtered_paths <- find_downstream_paths(
  start_nodes = c("Behaviorally plastic", "Physiologically plastic"),
  all_paths = filtered_coefficients_plasticity_mapped
)

# Extract unique display labels for nodes involved in these paths
unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

# Initialize the DOT script
dot_script <- "digraph SEM_PathDiagram { 
  rankdir=LR;
  node [shape=rectangle, style=filled, fillcolor=lightblue];
  
"

# Add nodes with shape customization and background color
for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    # Use rectangle shape with dark gray background for baseline nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    # Use default rectangle shape with white background for other nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

# Add edges with labels, colors, and penwidth
for(i in 1:nrow(filtered_paths)){
  predictor <- filtered_paths$Predictor_Display_Label[i]
  response <- filtered_paths$Response_Display_Label[i]
  estimate <- filtered_paths$Estimate[i]
  color <- get_edge_color(estimate)
  penwidth <- filtered_paths$Penwidth[i]
  label <- sprintf("%.2f", estimate)
  
  dot_script <- paste0(dot_script, sprintf("  \"%s\" -> \"%s\" [color=%s, label=\"%s\", penwidth=%.2f];\n",
                                           predictor, response, color, label, penwidth))
}

# Close the DOT script
dot_script <- paste0(dot_script, "}")

# Render the graph using grViz
grViz(dot_script)


### Pretty SEM Plots ----

## Herbivory SEM (Pretty) ----

# Define custom colors for treatments
treatment_colors <- c(
  "Herbivore" = "#B5753C",    # Brown
  "Vegetation" = "#9BA48C"     # Green
)

# Calculate penwidths (using the largest absolute effect as reference)
scale_penwidth <- function(value) {
  abs_value <- abs(value)
  min_width <- 1
  max_width <- 7
  max_effect <- 46.49  # Largest effect in the herbivory model
  
  scaled <- (abs_value / max_effect) * (max_width - min_width) + min_width
  return(scaled)
}

treatment_dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB;
  node [shape=rectangle, style=filled, fontname=\"Arial\"];
  splines=true;
  
  # Define node ranks
  { rank=source; \"Baseline goldenrod biomass\"; \"Herbivore\"; \"Baseline grass biomass\"; \"Vegetation\"; \"Baseline soil %N\" }
  { rank=same; \"Goldenrod biomass\"; \"Grass biomass\"; \"Soil %N\" }
  { rank=same; \"Plant diversity\"; \"SIR\" }
  { rank=same; \"Grass %C\"; \"Goldenrod %N\"; \"Nitrogen mineralization\" }
  { rank=sink; \"Soil %C\" }
  
  # Treatment nodes
  \"Herbivore\" [fillcolor=\"#B5753C\", fontcolor=\"white\"];
  \"Vegetation\" [fillcolor=\"#9BA48C\", fontcolor=\"white\"];
  
  # Baseline nodes
  \"Baseline goldenrod biomass\" [fillcolor=\"#A8A8A8\"];
  \"Baseline grass biomass\" [fillcolor=\"#A8A8A8\"];
  \"Baseline soil %N\" [fillcolor=\"#A8A8A8\"];
  
  # Response nodes
  \"Grass biomass\" [fillcolor=\"#F0F0F0\"];
  \"Goldenrod biomass\" [fillcolor=\"#F0F0F0\"];
  \"Soil %N\" [fillcolor=\"#F0F0F0\"];
  \"Soil %C\" [fillcolor=\"#F0F0F0\"];
  \"Plant diversity\" [fillcolor=\"#F0F0F0\"];
  \"Grass %C\" [fillcolor=\"#F0F0F0\"];
  \"Goldenrod %N\" [fillcolor=\"#F0F0F0\"];
  \"SIR\" [fillcolor=\"#F0F0F0\"];
  \"Nitrogen mineralization\" [fillcolor=\"#F0F0F0\"];
  
  # Treatment effects - Herbivore (Brown)
  \"Herbivore\" -> \"Goldenrod biomass\" [color=\"#B5753C\", label=\" 43.17\", penwidth=6.51];
  \"Herbivore\" -> \"Grass biomass\" [color=\"#B5753C\", label=\" 46.49\", penwidth=7.00];
  \"Herbivore\" -> \"Nitrogen mineralization\" [color=\"#B5753C\", label=\" 0.52\", penwidth=1.08];
  
  # Treatment effects - Vegetation (Green)
  \"Vegetation\" -> \"Grass biomass\" [color=\"#9BA48C\", label=\" 39.64\", penwidth=6.02];
  \"Vegetation\" -> \"Soil %N\" [color=\"#9BA48C\", label=\" 0.25\", penwidth=1.04];
  \"Vegetation\" -> \"Nitrogen mineralization\" [color=\"#9BA48C\", label=\" 13.39\", penwidth=2.73];
  
  # Baseline effects
  \"Baseline soil %N\" -> \"Soil %N\" [color=\"black\", label=\" 2.30\", penwidth=1.35];
  \"Baseline grass biomass\" -> \"Grass biomass\" [color=\"black\", label=\" 7.54\", penwidth=1.97];
  \"Baseline goldenrod biomass\" -> \"Goldenrod biomass\" [color=\"black\", label=\" 10.55\", penwidth=2.36];
  
  # Downstream effects (brown color)
  \"Grass biomass\" -> \"Plant diversity\" [color=\"#A18D6B\", label=\" -0.01\", penwidth=1.00, style=\"dashed\"];
  \"Grass biomass\" -> \"Grass %C\" [color=\"#A18D6B\", label=\" 0.08\", penwidth=1.01];
  \"Goldenrod biomass\" -> \"Plant diversity\" [color=\"#A18D6B\", label=\" -0.09\", penwidth=1.01, style=\"dashed\"];
  \"Goldenrod biomass\" -> \"Goldenrod %N\" [color=\"#A18D6B\", label=\" 0.01\", penwidth=1.00];
  \"Soil %N\" -> \"SIR\" [color=\"#A18D6B\", label=\" 3.23\", penwidth=1.46];
  \"SIR\" -> \"Soil %C\" [color=\"#A18D6B\", label=\" 0.64\", penwidth=1.10];
}
"

# Render the plot
grViz(treatment_dot_script)

## Plasticity SEM (Pretty) ----

# Define custom colors for plasticity types
plasticity_colors <- c(
  "Reactors" = "#663366",    # Muted Purple
  "Resistors" = "#D4AF37"    # Muted Gold
)

# Calculate penwidths (using the largest absolute effect of 10.50 as reference)
# Function to scale penwidth between 1 and 7
scale_penwidth <- function(value) {
  abs_value <- abs(value)
  min_width <- 1
  max_width <- 7
  max_effect <- 10.50
  
  scaled <- (abs_value / max_effect) * (max_width - min_width) + min_width
  return(scaled)
}

plasticity_dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB;
  node [shape=rectangle, style=filled, fontname=\"Arial\"];
  splines=true;
  
  # Define node ranks
  { rank=source; \"Baseline goldenrod biomass\"; \"Reactors\"; \"Resistors\"; \"Baseline soil %C\" }
  { rank=same; \"Goldenrod biomass\"; \"Grass biomass\"; \"Soil %N\" }
  { rank=same; \"Plant diversity\"; \"SIR\" }
  { rank=same; \"Grass %C\"; \"Goldenrod %N\"; \"Nitrogen mineralization\" }
  { rank=sink; \"Soil %C\" }
  
  # Treatment nodes
  \"Reactors\" [fillcolor=\"#663366\", fontcolor=\"white\"];
  \"Resistors\" [fillcolor=\"#D4AF37\", fontcolor=\"black\"];
  
  # Baseline nodes
  \"Baseline goldenrod biomass\" [fillcolor=\"#A8A8A8\"];
  \"Baseline soil %C\" [fillcolor=\"#A8A8A8\"];
  
  # Response nodes
  \"Grass biomass\" [fillcolor=\"#F0F0F0\"];
  \"Goldenrod biomass\" [fillcolor=\"#F0F0F0\"];
  \"Soil %N\" [fillcolor=\"#F0F0F0\"];
  \"Soil %C\" [fillcolor=\"#F0F0F0\"];
  \"Plant diversity\" [fillcolor=\"#F0F0F0\"];
  \"Grass %C\" [fillcolor=\"#F0F0F0\"];
  \"Goldenrod %N\" [fillcolor=\"#F0F0F0\"];
  \"SIR\" [fillcolor=\"#F0F0F0\"];
  \"Nitrogen mineralization\" [fillcolor=\"#F0F0F0\"];
  
  # Treatment effects - Reactors
  \"Reactors\" -> \"Goldenrod biomass\" [color=\"#663366\", label=\" -7.36\", penwidth=5.91, style=\"dashed\"];
  \"Reactors\" -> \"Grass biomass\" [color=\"#663366\", label=\" 6.16\", penwidth=5.12];
  \"Reactors\" -> \"Soil %N\" [color=\"#663366\", label=\" -0.02\", penwidth=1.02, style=\"dashed\"];
  \"Reactors\" -> \"Soil %C\" [color=\"#663366\", label=\" -0.16\", penwidth=1.11, style=\"dashed\"];
  \"Reactors\" -> \"Nitrogen mineralization\" [color=\"#663366\", label=\" 1.32\", penwidth=1.88];
  
  # Treatment effects - Resistors
  \"Resistors\" -> \"Goldenrod biomass\" [color=\"#D4AF37\", label=\" -4.68\", penwidth=4.13, style=\"dashed\"];
  \"Resistors\" -> \"Grass biomass\" [color=\"#D4AF37\", label=\" 7.54\", penwidth=6.02];
  \"Resistors\" -> \"Soil %N\" [color=\"#D4AF37\", label=\" -0.02\", penwidth=1.02, style=\"dashed\"];
  \"Resistors\" -> \"Soil %C\" [color=\"#D4AF37\", label=\" -0.19\", penwidth=1.13, style=\"dashed\"];
  \"Resistors\" -> \"Nitrogen mineralization\" [color=\"#D4AF37\", label=\" -3.07\", penwidth=2.95, style=\"dashed\"];
  
  # Baseline effects
  \"Baseline goldenrod biomass\" -> \"Goldenrod biomass\" [color=\"black\", label=\" 10.50\", penwidth=7.00];
  \"Baseline soil %C\" -> \"Grass %C\" [color=\"black\", label=\" 1.09\", penwidth=1.73];
  
  # Downstream effects
  \"Goldenrod biomass\" -> \"Goldenrod %N\" [color=\"#9F8151\", label=\" 0.01\", penwidth=1.01];
  \"Goldenrod biomass\" -> \"Plant diversity\" [color=\"#9F8151\", label=\" 0.00\", penwidth=1.00, style=\"dashed\"];
  \"Grass biomass\" -> \"Plant diversity\" [color=\"#9F8151\", label=\" -0.01\", penwidth=1.01, style=\"dashed\"];
  \"Grass biomass\" -> \"Grass %C\" [color=\"#9F8151\", label=\" 0.08\", penwidth=1.05];
  \"Soil %N\" -> \"SIR\" [color=\"#9F8151\", label=\" 3.23\", penwidth=3.06];
  \"SIR\" -> \"Soil %C\" [color=\"#9F8151\", label=\" 0.64\", penwidth=1.43];
}
"

# Render the plot
grViz(plasticity_dot_script)
