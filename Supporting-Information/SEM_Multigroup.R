# Supplement: Multigroup Piecewise SEM Analysis

library(lme4)
library(piecewiseSEM)
library(dplyr)

# Assumes final_data_year is already loaded and preprocessed as in SEM_Composite.R

# Split data by Treatment_PopType
data_veg <- final_data_year %>% filter(Treatment_PopType == "Vegetation")
data_phys <- final_data_year %>% filter(Treatment_PopType == "Reactor")
data_behav <- final_data_year %>% filter(Treatment_PopType == "Resistor")

# --- Group size and missingness check ---
cat("Vegetation group size:", nrow(data_veg), "\n")
cat("Reactor group size:", nrow(data_phys), "\n")
cat("Resistor group size:", nrow(data_behav), "\n\n")

cat("Missingness summary for Vegetation group:\n")
print(sapply(data_veg, function(x) sum(is.na(x))))
cat("\nMissingness summary for Reactor group:\n")
print(sapply(data_phys, function(x) sum(is.na(x))))
cat("\nMissingness summary for Resistor group:\n")
print(sapply(data_behav, function(x) sum(is.na(x))))
cat("\n")

# no missing data for groups 


# Helper function to fit models, dropping random effect if singular
fit_mixed_or_fixed <- function(formula_mixed, formula_fixed, data) {
  m <- tryCatch(lmer(formula_mixed, data = data), error = function(e) NULL)
  if (!is.null(m) && !lme4::isSingular(m, tol = 1e-4)) {
    return(m)
  } else {
    m_fixed <- tryCatch(lm(formula_fixed, data = data), error = function(e) NULL)
    if (!is.null(m_fixed)) {
      message("Fitting as fixed-effect model due to singular fit or error.")
      return(m_fixed)
    } else {
      message("Model could not be fit (even as fixed-effect).")
      return(NULL)
    }
  }
}

group_models <- function(data) {
  list(
    SORU_model = fit_mixed_or_fixed(
      SORU_Biomass_2023 ~ PercentN_SOIL_2021 + POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021 + (1 | Site),
      SORU_Biomass_2023 ~ PercentN_SOIL_2021 + POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021,
      data
    ),
    POPRC_model = fit_mixed_or_fixed(
      POPRC_Biomass_2023 ~ POPRC_Biomass_2021 + average_SOM_2021 + (1 | Site),
      POPRC_Biomass_2023 ~ POPRC_Biomass_2021 + average_SOM_2021,
      data
    ),
    SORUN_model = fit_mixed_or_fixed(
      PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + PercentN_SORU_2021 + (1 | Site),
      PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + PercentN_SORU_2021,
      data
    ),
    SORUC_model = fit_mixed_or_fixed(
      PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + PercentC_SORU_2021 + (1 | Site),
      PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + PercentC_SORU_2021,
      data
    ),
    POPRCC_model = fit_mixed_or_fixed(
      PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + PlantDiversity_2021 + Overall_mineralization_rate_2021 + POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021 + (1 | Site),
      PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + PlantDiversity_2021 + Overall_mineralization_rate_2021 + POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021,
      data
    ),
    POPRCN_model = fit_mixed_or_fixed(
      PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + Overall_mineralization_rate_2021 + PercentN_POPRC_2021 + (1 | Site),
      PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + Overall_mineralization_rate_2021 + PercentN_POPRC_2021,
      data
    ),
    SoilC_model = fit_mixed_or_fixed(
      PercentC_SOIL_2023 ~ SORU_Biomass_2023 + CO2CperHourperg_2023 + CO2CperHourperg_2021 + PercentC_POPRC_2021 + average_SOM_2021 + PercentC_SOIL_2021 + (1 | Site),
      PercentC_SOIL_2023 ~ SORU_Biomass_2023 + CO2CperHourperg_2023 + CO2CperHourperg_2021 + PercentC_POPRC_2021 + average_SOM_2021 + PercentC_SOIL_2021,
      data
    ),
    SoilN_model = fit_mixed_or_fixed(
      PercentN_SOIL_2023 ~ PercentN_LITTER_2023 + PercentN_SOIL_2021 + average_SOM_2021 + CO2CperHourperg_2021 + (1 | Site),
      PercentN_SOIL_2023 ~ PercentN_LITTER_2023 + PercentN_SOIL_2021 + average_SOM_2021 + CO2CperHourperg_2021,
      data
    ),
    SIR_model = fit_mixed_or_fixed(
      CO2CperHourperg_2023 ~ PercentN_SOIL_2023 + average_SOM_2021 + PlantDiversity_2021 + PercentN_SOIL_2021 + PercentN_POPRC_2021 + POPRC_Biomass_2021 + CO2CperHourperg_2021 + PercentN_POPRC_2021 + (1 | Site),
      CO2CperHourperg_2023 ~ PercentN_SOIL_2023 + average_SOM_2021 + PlantDiversity_2021 + PercentN_SOIL_2021 + PercentN_POPRC_2021 + POPRC_Biomass_2021 + CO2CperHourperg_2021 + PercentN_POPRC_2021,
      data
    ),
    PlantDiversity_model = fit_mixed_or_fixed(
      PlantDiversity_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + POPRC_Biomass_2023 + PercentN_SORU_2021 + average_SOM_2021 + PlantDiversity_2021 + Overall_mineralization_rate_2021 + PercentN_LITTER_2021 + (1 | Site),
      PlantDiversity_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + POPRC_Biomass_2023 + PercentN_SORU_2021 + average_SOM_2021 + PlantDiversity_2021 + Overall_mineralization_rate_2021 + PercentN_LITTER_2021,
      data
    ),
    LitterN_model = fit_mixed_or_fixed(
      PercentN_LITTER_2023 ~ PercentN_POPRC_2023 + PercentN_SORU_2023 + PercentN_LITTER_2021 + (1 | Site),
      PercentN_LITTER_2023 ~ PercentN_POPRC_2023 + PercentN_SORU_2023 + PercentN_LITTER_2021,
      data
    ),
    Nmin_model = fit_mixed_or_fixed(
      Overall_mineralization_rate_2023 ~ PercentN_SOIL_2023 + CO2CperHourperg_2023 + PercentN_POPRC_2021 + Overall_mineralization_rate_2021 + (1 | Site),
      Overall_mineralization_rate_2023 ~ PercentN_SOIL_2023 + CO2CperHourperg_2023 + PercentN_POPRC_2021 + Overall_mineralization_rate_2021,
      data
    )
  )
}

# Fit models
models_veg <- group_models(data_veg)
models_phys <- group_models(data_phys)
models_behav <- group_models(data_behav)

# Build SEMs for each group
sem_veg <- psem(
  models_veg$SORU_model,
  models_veg$POPRC_model,
  models_veg$SORUN_model,
  models_veg$SORUC_model,
  models_veg$POPRCC_model,
  models_veg$POPRCN_model,
  models_veg$SoilC_model,
  models_veg$SoilN_model,
  models_veg$SIR_model,
  models_veg$PlantDiversity_model,
  models_veg$LitterN_model,
  models_veg$Nmin_model,
  data = data_veg
)

sem_phys <- psem(
  models_phys$SORU_model,
  models_phys$POPRC_model,
  models_phys$SORUN_model,
  models_phys$SORUC_model,
  models_phys$POPRCC_model,
  models_phys$POPRCN_model,
  models_phys$SoilC_model,
  models_phys$SoilN_model,
  models_phys$SIR_model,
  models_phys$PlantDiversity_model,
  models_phys$LitterN_model,
  models_phys$Nmin_model,
  data = data_phys
)

sem_behav <- psem(
  models_behav$SORU_model,
  models_behav$POPRC_model,
  models_behav$SORUN_model,
  models_behav$SORUC_model,
  models_behav$POPRCC_model,
  models_behav$POPRCN_model,
  models_behav$SoilC_model,
  models_behav$SoilN_model,
  models_behav$SIR_model,
  models_behav$PlantDiversity_model,
  models_behav$LitterN_model,
  models_behav$Nmin_model,
  data = data_behav
)

sink("SEM_Summary_Vegetation.txt")
print(summary(sem_veg))
sink()

sink("SEM_Summary_Herbivore_Physiological.txt")
print(summary(sem_phys))
sink()

sink("SEM_Summary_Herbivore_Behavioral.txt")
print(summary(sem_behav))
sink() 
