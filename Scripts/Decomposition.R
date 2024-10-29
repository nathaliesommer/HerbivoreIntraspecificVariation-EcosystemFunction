library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(forcats)

### Load and clean data ----

# Read in data sets
decomp_2021 <- read.csv("Data/Decomposition/2021_Mesocosm-Decomp.csv")
decomp_2023 <- read.csv("Data/Decomposition/2023_Mesocosm-Decomp.csv")
field_litter <- read.csv("Data/Decomposition/2021_Field-Decomp.csv")

# Update columns with treatment types
decomp_2021 <- decomp_2021 %>%
  select(-Notes) %>%
  mutate(Site = case_when(
    startsWith(as.character(CageID), "FNH") ~ "FN",
    startsWith(as.character(CageID), "FNS") ~ "YF",
    startsWith(as.character(CageID), "FNN") ~ "UP",
    startsWith(as.character(CageID), "MCH") ~ "MC",
    startsWith(as.character(CageID), "MCN") ~ "UP",
    startsWith(as.character(CageID), "MCS") ~ "YF",
    startsWith(as.character(CageID), "YFH") ~ "YF",
    startsWith(as.character(CageID), "YFN") ~ "UP",
    startsWith(as.character(CageID), "SCN") ~ "UP",
    startsWith(as.character(CageID), "SCH") ~ "SC",
    startsWith(as.character(CageID), "SCS") ~ "YF",
    startsWith(as.character(CageID), "UPH") ~ "UP",
    startsWith(as.character(CageID), "UPS") ~ "YF",
    startsWith(as.character(CageID), "DCH") ~ "DC",
    startsWith(as.character(CageID), "DCN") ~ "UP",
    startsWith(as.character(CageID), "DCS") ~ "YF"
  )) %>%
  mutate(Latitude = case_when(
    Site %in% c("FN", "YF", "SC") ~ "South",
    Site %in% c("DC", "MC", "UP") ~ "North"
  )) %>%
  mutate(Collection_Weight = replace(Collection_Weight, Collection_Weight == 1222.000, 1.222)) %>%
  mutate(Trophic_Treatment = case_when(
    substr(CageID, 5, 5) == "V" ~ "Vegetation",
    substr(CageID, 5, 5) == "H" ~ "Herbivore",
    substr(CageID, 5, 5) == "P" ~ "Predator"
  )) %>%
  mutate(Transplant_Treatment = case_when(
    substr(CageID, 3, 3) == "H" ~ "Home",
    substr(CageID, 3, 3) == "N" ~ "North",
    substr(CageID, 3, 3) == "S" ~ "South"
  )) %>%
  mutate(Population = substr(CageID, 1, 2))

decomp_2023 <- decomp_2023 %>%
  select(-Notes) %>%
  mutate(Site = case_when(
    startsWith(as.character(CageID), "FNH") ~ "FN",
    startsWith(as.character(CageID), "FNS") ~ "YF",
    startsWith(as.character(CageID), "FNN") ~ "UP",
    startsWith(as.character(CageID), "MCH") ~ "MC",
    startsWith(as.character(CageID), "MCN") ~ "UP",
    startsWith(as.character(CageID), "MCS") ~ "YF",
    startsWith(as.character(CageID), "YFH") ~ "YF",
    startsWith(as.character(CageID), "YFN") ~ "UP",
    startsWith(as.character(CageID), "SCN") ~ "UP",
    startsWith(as.character(CageID), "SCH") ~ "SC",
    startsWith(as.character(CageID), "SCS") ~ "YF",
    startsWith(as.character(CageID), "UPH") ~ "UP",
    startsWith(as.character(CageID), "UPS") ~ "YF",
    startsWith(as.character(CageID), "DCH") ~ "DC",
    startsWith(as.character(CageID), "DCN") ~ "UP",
    startsWith(as.character(CageID), "DCS") ~ "YF"
  )) %>%
  mutate(Latitude = case_when(
    Site %in% c("FN", "YF", "SC") ~ "South",
    Site %in% c("DC", "MC", "UP") ~ "North"
  )) %>%
  mutate(Trophic_Treatment = case_when(
    substr(CageID, 5, 5) == "V" ~ "Vegetation",
    substr(CageID, 5, 5) == "H" ~ "Herbivore",
    substr(CageID, 5, 5) == "P" ~ "Predator"
  )) %>%
  mutate(Transplant_Treatment = case_when(
    substr(CageID, 3, 3) == "H" ~ "Home",
    substr(CageID, 3, 3) == "N" ~ "North",
    substr(CageID, 3, 3) == "S" ~ "South"
  )) %>%
  mutate(Population = substr(CageID, 1, 2))

### Calculate decomposition rates ----

# Litter-Standard Pairs
field_litter <- field_litter %>% 
  mutate(MassLoss = Set_Weight - Collection_Weight) %>%  
  rename(Site = Field) %>% 
  mutate(SetTime = as.numeric(difftime(mdy(Collection_Date), mdy(Set_Date), units = "days"))) %>% 
  mutate(DecompRate = (MassLoss / abs(SetTime)) * 30)

# Diagnostic plot: Comparison between standard and real litter
field_litter %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  ggplot(aes(y=DecompRate, x=Bag_Type)) +
  geom_boxplot(alpha=0.3, fill = "light gray") +  
  geom_jitter(width = .2, fill = "black", size = 1.1, alpha = 0.3) +
  theme_classic() +
  coord_flip() +
  labs(x = "",y = "Decomposition Rate (g/month)") + 
  facet_wrap(~Site)

# Calculate scaling factor for each site
scaling_factors <- field_litter %>%
  group_by(Site) %>%
  summarise(
    Real_Mean = mean(DecompRate[Bag_Type == "Litter"], na.rm = TRUE),
    Standard_Mean = mean(DecompRate[Bag_Type == "Standard"], na.rm = TRUE),
    Scaling_Factor = Real_Mean / Standard_Mean
  )

 # 2021 Cage-Level
litter_calc_2021 <- decomp_2021 %>%
  mutate(MassLoss = Set_Weight - Collection_Weight) %>%
  mutate(SetTime = as.numeric(difftime(mdy(Collection_Date), mdy(Set_Date), units = "days"))) %>%
  mutate(DecompRate = (MassLoss / abs(SetTime)) * 30) %>%
  mutate(Year = year(mdy(Collection_Date)))

# 2023 Cage-Level

# Calculate MassLoss for non-NA rows (2023 only)
decomp_2023 <- decomp_2023 %>%
  mutate(MassLoss = ifelse(!is.na(Set_Weight) & !is.na(Collection_Weight), Set_Weight - Collection_Weight, NA))

# Calculate average MassLoss for each combination of Site, Trophic_Treatment, and Population
average_massloss <- decomp_2023 %>%
  group_by(Site, Trophic_Treatment, Population) %>%
  summarise(Average_MassLoss = mean(MassLoss, na.rm = TRUE), .groups = 'drop')

# Interpolate NA MassLoss values
decomp_2023 <- decomp_2023 %>%
  left_join(average_massloss, by = c("Site", "Trophic_Treatment", "Population")) %>%
  mutate(MassLoss = ifelse(is.na(MassLoss), Average_MassLoss, MassLoss)) %>%
  select(-Average_MassLoss)

# Calculate decomposition rates
litter_calc_2023 <- decomp_2023 %>%
  mutate(SetTime = as.numeric(difftime(mdy(Collection_Date), mdy(Set_Date), units = "days"))) %>%
  mutate(DecompRate = (MassLoss / abs(SetTime)) * 30) %>%
  mutate(Year = year(mdy(Collection_Date)))

# Average replicates within each year 
decomp_2021_avg <- litter_calc_2021 %>%
  group_by(Site, CageID, Population, Trophic_Treatment) %>%
  summarise(Unscaled_MassLoss_2021 = mean(MassLoss, na.rm = TRUE), .groups = 'drop')

decomp_2023_avg <- decomp_2023 %>%
  group_by(Site, CageID, Population, Trophic_Treatment) %>%
  summarise(Unscaled_MassLoss_2023 = mean(MassLoss, na.rm = TRUE), .groups = 'drop')

# Rename columns before combining
decomp_2021_avg <- decomp_2021_avg %>%
  rename(MassLoss = Unscaled_MassLoss_2021) %>%
  mutate(Year = 2021)

decomp_2023_avg <- decomp_2023_avg %>%
  rename(MassLoss = Unscaled_MassLoss_2023) %>%
  mutate(Year = 2023)

# Combine the datasets
combined_decomp <- bind_rows(decomp_2021_avg, decomp_2023_avg)


