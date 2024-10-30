### N-mineralization ----
## Nitrate = (mg N-NO3/mL)
## Ammonium = (mg N-NH4/mL)
## Bulk density = g/cm^3

## Notes
# For those that remain below zero, still use these values as-is 
# in analyses of both nitrogen concentrations and in calculating 
# mineralization and nitrification rates

## Nitrification was calculated as:
# NO3 at 28 days minus initial NO3 at time zero; 
# N mineralization equals NH4 + NO3 at 28 days minus NH4 + NO3 at time zero (Goodale and Aber 2001).
# This was multiplied by bulk density and the appropriate conversions to obtain results on a g/cm^3 basis Pastor et al 1987: 


library(dplyr)

# Read-in data
N_min_2021 <- read.csv("Data/N-min/2021_N-min.csv")
N_min_2023 <- read.csv("Data/N-min/2023_N-min.csv") # NOTE: 2023 N-min data was blank-adjusted by the processing lab
BulkDensity <- read.csv("Data/BulkDensity/2021_BulkDensity.csv")

# Bulk Density calculation
BulkDensity <- BulkDensity %>% 
  mutate(SoilMass = Weight_SoilTray-Weight_Tray) %>% 
  mutate(SoilMassAdj = SoilMass - (Materials_volume_mL_end - Materials_volume_mL_start)) %>% 
  mutate(SoilBulkDensity = (SoilMassAdj / 308.89)) # units = g/cm^3, where 308.89 is the cm^3 of the sampler

# Add columns for treatment

N_min_2021 <- N_min_2021 %>%
  mutate(Site = case_when(
    startsWith(as.character(Sample.ID), "FNH") ~ "FN",
    startsWith(as.character(Sample.ID), "FNS") ~ "YF",
    startsWith(as.character(Sample.ID), "FNN") ~ "UP",
    startsWith(as.character(Sample.ID), "MCH") ~ "MC",
    startsWith(as.character(Sample.ID), "MCN") ~ "UP",
    startsWith(as.character(Sample.ID), "MCS") ~ "YF",
    startsWith(as.character(Sample.ID), "YFH") ~ "YF",
    startsWith(as.character(Sample.ID), "YFN") ~ "UP",
    startsWith(as.character(Sample.ID), "SCN") ~ "UP",
    startsWith(as.character(Sample.ID), "SCH") ~ "SC",
    startsWith(as.character(Sample.ID), "SCS") ~ "YF",
    startsWith(as.character(Sample.ID), "UPH") ~ "UP",
    startsWith(as.character(Sample.ID), "UPS") ~ "YF",
    startsWith(as.character(Sample.ID), "DCH") ~ "DC",
    startsWith(as.character(Sample.ID), "DCN") ~ "UP",
    startsWith(as.character(Sample.ID), "DCS") ~ "YF"
  )) %>%
  mutate(Trophic_Treatment = case_when(
    substr(Sample.ID, 5, 5) == "V" ~ "Vegetation",
    substr(Sample.ID, 5, 5) == "H" ~ "Herbivore",
    substr(Sample.ID, 5, 5) == "P" ~ "Predator"
  )) %>%
  mutate(Transplant_Treatment = case_when(
    substr(Sample.ID, 3, 3) == "H" ~ "Home",
    substr(Sample.ID, 3, 3) == "N" ~ "North",
    substr(Sample.ID, 3, 3) == "S" ~ "South"
  )) %>%
  mutate(Population = substr(Sample.ID, 1, 2))


N_min_2023 <- N_min_2023 %>%
  mutate(Site = case_when(
    startsWith(as.character(Sample.ID), "FNH") ~ "FN",
    startsWith(as.character(Sample.ID), "FNS") ~ "YF",
    startsWith(as.character(Sample.ID), "FNN") ~ "UP",
    startsWith(as.character(Sample.ID), "MCH") ~ "MC",
    startsWith(as.character(Sample.ID), "MCN") ~ "UP",
    startsWith(as.character(Sample.ID), "MCS") ~ "YF",
    startsWith(as.character(Sample.ID), "YFH") ~ "YF",
    startsWith(as.character(Sample.ID), "YFN") ~ "UP",
    startsWith(as.character(Sample.ID), "SCN") ~ "UP",
    startsWith(as.character(Sample.ID), "SCH") ~ "SC",
    startsWith(as.character(Sample.ID), "SCS") ~ "YF",
    startsWith(as.character(Sample.ID), "UPH") ~ "UP",
    startsWith(as.character(Sample.ID), "UPS") ~ "YF",
    startsWith(as.character(Sample.ID), "DCH") ~ "DC",
    startsWith(as.character(Sample.ID), "DCN") ~ "UP",
    startsWith(as.character(Sample.ID), "DCS") ~ "YF"
  )) %>%
  mutate(Trophic_Treatment = case_when(
    substr(Sample.ID, 5, 5) == "V" ~ "Vegetation",
    substr(Sample.ID, 5, 5) == "H" ~ "Herbivore",
    substr(Sample.ID, 5, 5) == "P" ~ "Predator"
  )) %>%
  mutate(Transplant_Treatment = case_when(
    substr(Sample.ID, 3, 3) == "H" ~ "Home",
    substr(Sample.ID, 3, 3) == "N" ~ "North",
    substr(Sample.ID, 3, 3) == "S" ~ "South"
  )) %>%
  mutate(Population = substr(Sample.ID, 1, 2))


# blank adjustment for 2021 data
N_min_2021 <- N_min_2021 %>%
  group_by(Batch) %>%
  mutate(
    N_NO3_blank_avg = mean(`N.NO3..mg.per.mL.`[Sample.ID == "Blank"]),
    N_NH4_blank_avg = mean(`N.NH4..mg.per.mL.`[Sample.ID == "Blank"])
  ) %>%
  mutate(
    `N-NO3 (mg per mL)` = `N.NO3..mg.per.mL.` - N_NO3_blank_avg,
    `N-NH4 (mg per mL)` = `N.NH4..mg.per.mL.` - N_NH4_blank_avg
  ) %>%
  ungroup() %>%
  select(-N_NO3_blank_avg, -N_NH4_blank_avg)

# Remove blanks from 2023 and 2021 data 
N_min_2023 <- N_min_2023 %>%
  filter(Sample.ID != "Blank")

N_min_2021 <- N_min_2021 %>% 
  filter(Sample.ID != "Blank")

# Bulk Density, unit conversion 
BulkDensity_avg <- BulkDensity %>%
  group_by(Site) %>%
  summarize(SoilBulkDensity_avg = mean(SoilBulkDensity))


calculate_rate <- function(data, year, soil_depth = 10) {
  data %>%
    left_join(BulkDensity_avg, by = "Site") %>%
    group_by(Sample.ID) %>%
    # Filter for the first and last day
    filter(Day == min(Day) | Day == max(Day)) %>%
    arrange(Day) %>%
    # Calculate differences between the last and first day
    mutate(
      NH4_diff = `N.NH4..mg.per.mL.`[2] - `N.NH4..mg.per.mL.`[1],
      NO3_diff = `N.NO3..mg.per.mL.`[2] - `N.NO3..mg.per.mL.`[1],
      # Calculate rates, scaled to per month and converted to mg/cm³
      NH4_rate = (NH4_diff * 1) * SoilBulkDensity_avg * soil_depth, # mg/cm³ per month
      NO3_rate = (NO3_diff * 1) * SoilBulkDensity_avg * soil_depth, # mg/cm³ per month
      Overall_rate = NH4_rate + NO3_rate
    ) %>%
    ungroup()
}


# Calculate rates and add the Year column
N_min_2021 <- calculate_rate(N_min_2021, year = 2021) %>%
  mutate(Year = 2021)

N_min_2023 <- calculate_rate(N_min_2023, year = 2023) %>%
  mutate(Year = 2023)

# Clean up the data frames after calculating rates
# Drop columns that exist only in some data frames
clean_data <- function(data, year_column_exists = TRUE) {
data <- data %>%
  select(-Day, -`N.NO3..mg.per.mL.`, -`N.NH4..mg.per.mL.`, 
         -`N.NO3..mg.per.mL.`, -`N.NH4..mg.per.mL.`, 
         -SoilBulkDensity_avg, -NH4_diff, -NO3_diff)

# Conditionally drop the Batch column if it exists
if ("Batch" %in% colnames(data)) {
  data <- data %>%
    select(-Batch)
}

# Ensure distinct rows based on Sample.ID and Year if Year column exists
if (year_column_exists) {
  data <- data %>%
    distinct(Sample.ID, Year, .keep_all = TRUE)
}

return(data)
}

N_min_2021_clean <- clean_data(N_min_2021)
N_min_2023_clean <- clean_data(N_min_2023)

# Merge the cleaned data frames
N_min_full_data <- bind_rows(N_min_2021_clean, N_min_2023_clean) %>%
  select(Sample.ID, Year, Transplant_Treatment, Trophic_Treatment, 
         Population, Site, NH4_rate, NO3_rate, Overall_rate) %>%
  rename(
    `Ammonium rate` = NH4_rate,
    `Nitrate rate` = NO3_rate,
    `Overall mineralization rate` = Overall_rate
  )

## Units are (mg N/cm² per month)





### Permute UPH_H7 and YFN_H2 2023

# Create new rows with specified information
new_rows <- tibble(
  Sample.ID = c("UPH_H7", "YFN_H2"),
  Year = c(2023, 2023),
  Transplant_Treatment = c("Home", "North"),
  Trophic_Treatment = c("Herbivore", "Herbivore"),
  Population = c("UP", "YF"),
  Site = c("UP", "UP")
)

# Function to calculate average rates based on matching conditions
average_rates <- function(data, new_row) {
  matching_rows <- data %>%
    filter(
      Year == new_row$Year,
      Population == new_row$Population,
      Site == new_row$Site,
      Trophic_Treatment == new_row$Trophic_Treatment
    )
  
  if (nrow(matching_rows) > 0) {
    Ammonium_rate <- mean(matching_rows$`Ammonium rate`, na.rm = TRUE)
    Nitrate_rate <- mean(matching_rows$`Nitrate rate`, na.rm = TRUE)
    Overall_mineralization_rate <- mean(matching_rows$`Overall mineralization rate`, na.rm = TRUE)
  } else {
    Ammonium_rate <- NA
    Nitrate_rate <- NA
    Overall_mineralization_rate <- NA
  }
  
  return(c(Ammonium_rate, Nitrate_rate, Overall_mineralization_rate))
}

# Apply the average calculation function to each new row
new_rows <- new_rows %>%
  rowwise() %>%
  mutate(
    rates = list(average_rates(N_min_full_data, as.list(cur_data()))),
    `Ammonium rate` = rates[[1]],
    `Nitrate rate` = rates[[2]],
    `Overall mineralization rate` = rates[[3]]
  ) %>%
  select(-rates) %>%
  ungroup()

# Combine the new rows with the existing data
N_min_full_data <- bind_rows(N_min_full_data, new_rows)

# View the updated data
print(N_min_full_data)








## GRAPHS ----
# Load required library
library(ggplot2)
library(tidyr)

# Filter the 2021 data
N_min_2021_clean <- N_min_full_data %>%
  filter(Year == 2021)

# Create a plot of overall rates by site
ggplot(N_min_2021_clean, aes(x = Site, y = `Overall mineralization rate`, fill = Site)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) + # Adds points for individual data
  theme_classic() +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() + 
  labs(
    title = "Total Nitrogen Mineralization Rate by Site",
    x = "Site",
    y = "Total Mineralization Rate (g N/cm² per month)"
  )

# Calculate the standard deviation of the overall rate by site
site_variability <- N_min_2021_clean %>%
  group_by(Site) %>%
  summarize(
    Within_Site_SD = sd(`Overall mineralization rate`, na.rm = TRUE),
    Sample_Size = n()
  )

# Calculate the overall standard deviation across all sites
overall_variability <- tibble(
  Site = "Overall",
  Within_Site_SD = sd(N_min_2021_clean$`Overall mineralization rate`, na.rm = TRUE),
  Sample_Size = nrow(N_min_2021_clean)
)

# Combine the site-specific variability with overall variability
variability_data <- bind_rows(site_variability, overall_variability)

# Add a column to distinguish between site-specific and overall variability
variability_data <- variability_data %>%
  mutate(Variability_Type = ifelse(Site == "Overall", "Overall Variability", "Within Site Variability"))

# Order the bars sequentially by SD
variability_data <- variability_data %>%
  arrange(desc(Within_Site_SD))

# Create the plot
ggplot(variability_data, aes(x = reorder(Site, -Within_Site_SD), y = Within_Site_SD, fill = Variability_Type)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_text(aes(label = Sample_Size), 
            position = position_dodge(width = 0.9), 
            vjust = -0.2, hjust = -0.2) +
  theme_classic() +
  labs(x = "Site", y = "Standard Deviation of Total Mineralization Rate (g N/cm² per month)") +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  ggtitle("Variability of Total Nitrogen Mineralization Within and Across Sites")


#### N-min differences -----
# Calculate the difference in overall nitrogen mineralization rate between 2023 and 2021
N_min_diff <- N_min_full_data %>%
  select(Sample.ID, Year, Population, Site, Transplant_Treatment, Trophic_Treatment, `Overall mineralization rate`) %>%
  spread(key = Year, value = `Overall mineralization rate`) %>%
  mutate(Diff_Overall_Rate = `2023` - `2021`) %>%
  select(Sample.ID, Population, Site, Transplant_Treatment, Trophic_Treatment, Diff_Overall_Rate) %>% 
  # Filter out rows where Site or Population is "MC"
  filter(Site != "MC", Population != "MC")

# Plot 1: Grouped by Trophic Treatment
ggplot(N_min_diff, aes(x = Population, y = Diff_Overall_Rate, fill = Trophic_Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +  # Align boxplots
  geom_jitter(aes(color = Trophic_Treatment), alpha = 0.5, position = position_dodge(width = 0.75)) +  # Align jitter points
  theme_classic() +
  labs(
    title = "Difference in Overall Nitrogen Mineralization Rate (Treatment - Baseline)",
    subtitle = "Grouped by Trophic Treatment",
    x = "Population",
    y = "Change in Overall Mineralization Rate (g N/cm² per month)",
    fill = "Trophic Treatment",
    color = "Trophic Treatment"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  coord_flip()

# Plot 2: Grouped by Transplant Treatment
ggplot(N_min_diff, aes(x = Population, y = Diff_Overall_Rate, fill = Transplant_Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.75)) +  # Align boxplots
  geom_jitter(aes(color = Transplant_Treatment), alpha = 0.5, position = position_dodge(width = 0.75)) +  # Align jitter points
  theme_classic() +
  labs(
    title = "Difference in Overall Nitrogen Mineralization Rate (Treatment - Baseline)",
    subtitle = "Grouped by Transplant Treatment",
    x = "Population",
    y = "Change in Overall Mineralization Rate (g N/cm² per month)",
    fill = "Transplant Treatment",
    color = "Transplant Treatment"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  coord_flip()


