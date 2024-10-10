# Load necessary libraries
library(tidyverse)

# Read the data files
gwc_2021 <- read_csv("Data/GravimetricWater/2021_GWC.csv")
gwc_2023 <- read_csv("Data/GravimetricWater/2023_GWC.csv")

# Define a function to perform the necessary mutations
mutate_gwc_data <- function(data) {
  data %>%
    mutate(
      Site = case_when(
        startsWith(as.character(Cage_ID), "FNH") ~ "FN",
        startsWith(as.character(Cage_ID), "FNS") ~ "YF",
        startsWith(as.character(Cage_ID), "FNN") ~ "UP",
        startsWith(as.character(Cage_ID), "MCH") ~ "MC",
        startsWith(as.character(Cage_ID), "MCN") ~ "UP",
        startsWith(as.character(Cage_ID), "MCS") ~ "YF",
        startsWith(as.character(Cage_ID), "YFH") ~ "YF",
        startsWith(as.character(Cage_ID), "YFN") ~ "UP",
        startsWith(as.character(Cage_ID), "SCN") ~ "UP",
        startsWith(as.character(Cage_ID), "SCH") ~ "SC",
        startsWith(as.character(Cage_ID), "SCS") ~ "YF",
        startsWith(as.character(Cage_ID), "UPH") ~ "UP",
        startsWith(as.character(Cage_ID), "UPS") ~ "YF",
        startsWith(as.character(Cage_ID), "DCH") ~ "DC",
        startsWith(as.character(Cage_ID), "DCN") ~ "UP",
        startsWith(as.character(Cage_ID), "DCS") ~ "YF"
      ),
      Trophic_Treatment = case_when(
        substr(Cage_ID, 5, 5) == "V" ~ "Vegetation",
        substr(Cage_ID, 5, 5) == "H" ~ "Herbivore",
        substr(Cage_ID, 5, 5) == "P" ~ "Predator"
      ),
      Transplant_Treatment = case_when(
        substr(Cage_ID, 3, 3) == "H" ~ "Home",
        substr(Cage_ID, 3, 3) == "N" ~ "North",
        substr(Cage_ID, 3, 3) == "S" ~ "South"
      ),
      Population = substr(Cage_ID, 1, 2)
    )
}

# Apply the mutations to both datasets
gwc_2021 <- mutate_gwc_data(gwc_2021)
gwc_2023 <- mutate_gwc_data(gwc_2023)

# Plot 1: Moisture across sites in 2021
ggplot(gwc_2021, aes(x = Site, y = `Moist (g H2O g soil-1)`)) +
  geom_boxplot() +
  labs(title = "Moisture across sites in 2021", x = "Site", y = "Moist (g H2O g soil -1)")

# Combine the datasets for change analysis
gwc_combined <- full_join(gwc_2021, gwc_2023, by = "Cage_ID", suffix = c("_2021", "_2023"))

# Calculate the change in moisture
gwc_combined <- gwc_combined %>%
  mutate(Moisture_Change = `Moist (g H2O g soil-1)_2023` - `Moist (g H2O g soil-1)_2021`)

# Drop missing or infinite values
gwc_combined <- gwc_combined %>%
  filter(!is.na(Moisture_Change) & is.finite(Moisture_Change))

# Plot 2: Change in moisture between 2021 and 2023 across Sites
ggplot(gwc_combined, aes(x = Site_2021, y = Moisture_Change)) +
  geom_boxplot() +
  labs(title = "Change in moisture between 2021 and 2023 across Sites", x = "Site", y = "Change in Moist (g H2O g soil -1)")

# Plot 3: Change in moisture between 2021 and 2023 across Trophic_Treatment
ggplot(gwc_combined, aes(x = Trophic_Treatment_2021, y = Moisture_Change)) +
  geom_boxplot() +
  labs(title = "Change in moisture between 2021 and 2023 across Trophic Treatment", x = "Trophic Treatment", y = "Change in Moist (g H2O g soil -1)")
