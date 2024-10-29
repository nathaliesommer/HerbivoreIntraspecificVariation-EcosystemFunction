# Vegetation Allometry

library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)

## Data Import and Cleaning ----

### Biomass ----
vegmass2021 <- read.csv("Data/VegBiomass/2021_Biomass-plots.csv") # biomass plots 2021
vegmass2023 <- read.csv("Data/VegBiomass/2023_Biomass-plots.csv") # biomass plots 2023

vegmass2023 <- vegmass2023 %>% # unit conversion
  mutate(
    Biomass_Bag_g = Biomass_Bag_kg * 1000,
    Bag_g = Bag_kg * 1000,
    Biomass_g = Biomass_kg * 1000
  )

vegmass2023 <- vegmass2023 %>%
  select(-Biomass_Bag_kg, -Bag_kg, -Biomass_kg) 

biomass_plots <- bind_rows(vegmass2021, vegmass2023) # bind 2021 and 2023 data sets

### Cover ----
diversity2021 <- read.csv("Data/VegCommunity/Vegetation.cover.2021.csv") # plant percent cover 2021
diversity2023 <- read.csv("Data/VegCommunity/Vegetation.cover.2023.csv") # plant percent cover 2023

diversitycomb <- bind_rows(diversity2023, diversity2021) # bind 2021 and 2023 data sets

cage_diversity <- diversitycomb %>%
  filter(!(Population %in% c("HF", "SP"))) %>% # remove extraneous sites
  select(-BARE) %>% # remove bare ground
  rowwise() %>%
  mutate(Cover_Sum = sum(c_across(-c(Cage.ID, Year, Population, Transplant, Treatment, Site, Replicate)))) %>% # sum plant cover 
  ungroup()

cage_diversity <- cage_diversity %>% # normalize cover to 100%
  mutate(across(
    .cols = -c(Cage.ID, Year, Population, Transplant, Treatment, Site, Replicate, Cover_Sum),
    .fns = ~ round((.x / Cover_Sum) * 100)
  ))

cage_diversity_long <- cage_diversity %>%
  pivot_longer(
    cols = -c(Cage.ID, Year, Population, Transplant, Treatment, Site, Replicate, Cover_Sum),
    names_to = "Species_ID",
    values_to = "Cover"
  ) %>% 
  select(-Cover_Sum)



functional_groups <- cage_diversity_long %>%
  group_by(Cage.ID, Year, Population, Treatment, Transplant, Site, Replicate) %>% 
  summarise(
    SORU = sum(Cover[Species_ID %in% c("SORU2", "SOCA6")], na.rm = TRUE),
    POPRC = sum(Cover[Species_ID == "POPRC"], na.rm = TRUE),
    MISC = sum(Cover[!(Species_ID %in% c("SORU2", "SOCA6", "POPRC"))], na.rm = TRUE)
  ) %>%
  pivot_longer(cols = SORU:MISC, names_to = "Species_ID", values_to = "Cover") %>%
  ungroup() %>% 
  drop_na()



## Allometry ----
# Use biomass plots to estimate the biomass in cages from % cover

# Calculate the biomass for each plant group within biomass_plots
biomass_plots <- biomass_plots %>%
  mutate(
    SORU_Biomass = (SORU / 100) * Biomass_g,
    POPRC_Biomass = (POPRC / 100) * Biomass_g,
    MISC_Biomass = (MISC / 100) * Biomass_g
  )

# Nest data by Site and Year
nested_data <- biomass_plots %>%
  group_by(Site, Year) %>%
  nest()

# Fit separate models and generate diagnostic plots for each group
model_fits <- nested_data %>%
  mutate(
    # Fit models for each group
    SORU_model = map(data, ~ lm(SORU_Biomass ~ SORU, data = .x)),
    POPRC_model = map(data, ~ lm(POPRC_Biomass ~ POPRC, data = .x)),
    MISC_model = map(data, ~ lm(MISC_Biomass ~ MISC, data = .x)),
    
    # Diagnostic plots for SORU
    diagnostic_plot_SORU = map2(data, Site, ~ ggplot(.x, aes(x = SORU, y = SORU_Biomass)) +
                                  geom_point() +
                                  geom_smooth(method = "lm", se = FALSE) +
                                  ggtitle(paste("Site:", .y, "SORU Check"))
                                ),
    # Diagnostic plots for POPRC
    diagnostic_plot_POPRC = map2(data, Site, ~ ggplot(.x, aes(x = POPRC, y = POPRC_Biomass)) +
                                   geom_point() +
                                   geom_smooth(method = "lm", se = FALSE) +
                                   ggtitle(paste("Site:", .y, "POPRC Check"))
                                ),
    
    # Diagnostic plots for MISC
    diagnostic_plot_MISC = map2(data, Site, ~ ggplot(.x, aes(x = MISC, y = MISC_Biomass)) +
                                  geom_point() +
                                  geom_smooth(method = "lm", se = FALSE) +
                                  ggtitle(paste("Site:", .y, "MISC Check"))
                                )
  )

# View diagnostic plots, check if transformations are needed
model_fits %>% 
  pull(diagnostic_plot_SORU) %>%
  walk(print)

model_fits %>% 
  pull(diagnostic_plot_POPRC) %>%
  walk(print)

model_fits %>% 
  pull(diagnostic_plot_MISC) %>%
  walk(print)



functional_groups_wide <- functional_groups %>%
  pivot_wider(names_from = Species_ID, values_from = Cover) %>%
  # Replace missing cover values (e.g., if some species are missing from some plots) with 0
  mutate(
    SORU = as.numeric(coalesce(SORU, 0)),
    POPRC = as.numeric(coalesce(POPRC, 0)),
    MISC = as.numeric(coalesce(MISC, 0))
  )

# Apply models row by row for each combination of Site and Year
functional_groups_wide <- functional_groups_wide %>%
  rowwise() %>%  # This ensures we work on each row individually
  mutate(
    # Predict SORU biomass
    SORU_Biomass = {
      model_idx <- which(model_fits$Site == Site & model_fits$Year == Year)
      predict(model_fits$SORU_model[[model_idx]], newdata = data.frame(SORU = SORU))
    },
    # Predict POPRC biomass
    POPRC_Biomass = {
      model_idx <- which(model_fits$Site == Site & model_fits$Year == Year)
      predict(model_fits$POPRC_model[[model_idx]], newdata = data.frame(POPRC = POPRC))
    },
    # Predict MISC biomass
    MISC_Biomass = {
      model_idx <- which(model_fits$Site == Site & model_fits$Year == Year)
      predict(model_fits$MISC_model[[model_idx]], newdata = data.frame(MISC = MISC))
    }
  ) %>%
  ungroup()

### Diversity Index Calculation ----

# Function to calculate Shannon-Weiner diversity index
calculate_shannon_weiner <- function(abundances) {
  proportions <- abundances / sum(abundances)
  proportions <- proportions[proportions > 0]  # Remove zero proportions
  -sum(proportions * log(proportions))
}

# Function to calculate diversity indices
calculate_diversity_indices <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      PlantRichness = sum(c_across(-c(Cage.ID, Year, Population, Transplant, Treatment, Site, Replicate)) > 0),
      PlantDiversity = calculate_shannon_weiner(c_across(-c(Cage.ID, Year, Population, Transplant, Treatment, Site, Replicate)))
    ) %>%
    ungroup() %>%
    select(Cage.ID, PlantRichness, PlantDiversity)
}

# Load your data
veg_2021 <- read.csv("Data/VegCommunity/Vegetation.cover.2021.csv")
veg_2023 <- read.csv("Data/VegCommunity/Vegetation.cover.2023.csv")

# Calculate diversity indices for each dataset
diversity_2021 <- calculate_diversity_indices(veg_2021) %>% 
  drop_na() # removes cages from data set not relevant to this study
diversity_2023 <- calculate_diversity_indices(veg_2023) %>% 
  drop_na() # removes cages from data set not relevant to this study

# Combine results and transform to long format
combined_diversity_long <- full_join(diversity_2021, diversity_2023, by = "Cage.ID", suffix = c("_2021", "_2023")) %>%
  pivot_longer(
    cols = starts_with("Plant"),
    names_to = c(".value", "Year"),
    names_sep = "_"
  )

