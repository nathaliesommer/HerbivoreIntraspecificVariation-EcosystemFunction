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
  select(Sample_ID, Population, Scaled_MassLoss, Year)

# SIR data
sir_data <- SIR_final_data %>%
  select(Sample_ID, Year, CO2CperHourperg)

# N-mineralization data
n_min_data <- N_min_full_data %>%
  rename(Sample_ID = Sample.ID) %>%
  select(Sample_ID, Year, Population, Site, `Ammonium rate`, `Nitrate rate`)

# Vegetation biomass data
veg_biomass_data <- functional_groups_wide %>%
  rename(Sample_ID = Cage.ID) %>%
  select(Sample_ID, Year, Population, Treatment, Transplant, Site, SORU_Biomass, POPRC_Biomass, MISC_Biomass)

# Plant diversity data
diversity_data <- combined_diversity_long %>%
  rename(Sample_ID = Cage.ID) %>%
  select(Sample_ID, Year, PlantRichness, PlantDiversity)

# Soil and plant CN data
cn_data <- CNdata %>%
  rename(Sample_ID = CageID) %>%
  select(Sample_ID, Year, Population, Site, PercentN, PercentC, SampleType) %>%
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
  select(
    Sample_ID,
    Year,
    Population = Population.x,  
    Site = Site.x,
    Treatment,
    Transplant,
    Scaled_MassLoss,
    CO2CperHourperg,
    `Ammonium rate`,
    `Nitrate rate`,
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
    )
  )

# Check the final aggregated data frame
head(final_data)



### Check assumptions ----

library(lme4)
library(ggplot2)
library(performance)
library(DHARMa)

#### Decomposition ----

Decomp_model <- lmer(Scaled_MassLoss ~ Year * Treatment * Population + (1|Site/Sample_ID), data = final_data)

# 1. DHARMa Diagnostics
simulation_output <- simulateResiduals(fittedModel = Decomp_model, n = 1000)
plot(simulation_output) # underdispersion
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # not zero-inflated

# 2. Multicollinearity Assessment
fixed_model <- lm(Scaled_MassLoss ~ Year * Treatment * Population, data = final_data)
vif_values <- vif(fixed_model, type = "predictor")
print(vif_values) 

vif_values <- vif(fixed_model)
print(vif_values) 

# Interaction effects are introducing collinearity; but the main effects are independent

# 3. Transformations
##### START HERE -----**

#### SIR ----

SIR_model <- lmer(CO2CperHourperg ~ Year * Treatment * Population + (1|Site/Sample_ID), data = final_data)

# 1. DHARMa Diagnostics
simulation_output <- simulateResiduals(fittedModel = SIR_model, n = 1000)
plot(simulation_output)
testDispersion(simulation_output) 
testZeroInflation(simulation_output) # not zero-inflated

# 2. Multicollinearity Assessment
fixed_model <- lm(CO2CperHourperg ~ Year * Treatment * Population, data = final_data)
vif_values <- vif(fixed_model, type = "predictor")
print(vif_values)

vif_values <- vif(fixed_model)
print(vif_values)
