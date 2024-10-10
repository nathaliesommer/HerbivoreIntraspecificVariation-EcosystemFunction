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

 # 2021 Cage-Level
litter_calc_2021 <- decomp_2021 %>%
  mutate(MassLoss = Set_Weight - Collection_Weight) %>%
  mutate(SetTime = as.numeric(difftime(mdy(Collection_Date), mdy(Set_Date), units = "days"))) %>%
  mutate(DecompRate = (MassLoss / abs(SetTime)) * 30) %>%
  mutate(Year = year(mdy(Collection_Date)))

# 2023 Cage-Level
litter_calc_2023 <- decomp_2023 %>%
  mutate(MassLoss = Set_Weight - Collection_Weight) %>%
  mutate(SetTime = as.numeric(difftime(mdy(Collection_Date), mdy(Set_Date), units = "days"))) %>%
  mutate(DecompRate = (MassLoss / abs(SetTime)) * 30) %>%
  mutate(Year = year(mdy(Collection_Date)))

### Population Effects: 2023 and 2021 comparison ----

# Calculate scaling factor for each site
scaling_factors <- field_litter %>%
  group_by(Site) %>%
  summarise(
    Real_Mean = mean(DecompRate[Bag_Type == "Litter"], na.rm = TRUE),
    Standard_Mean = mean(DecompRate[Bag_Type == "Standard"], na.rm = TRUE),
    Scaling_Factor = Real_Mean / Standard_Mean
  )

# Apply scaling factors to 2021 and 2023 data
decomp_2021_scaled <- litter_calc_2021 %>%
  left_join(scaling_factors, by = "Site") %>%
  mutate(Scaled_MassLoss = MassLoss * Scaling_Factor) %>%
  select(Site, CageID, Population, Trophic_Treatment, Transplant_Treatment, Rep, Scaled_MassLoss)

decomp_2023_scaled <- litter_calc_2023 %>%
  left_join(scaling_factors, by = "Site") %>%
  mutate(Scaled_MassLoss = MassLoss * Scaling_Factor) %>%
  select(Site, CageID, Population, Trophic_Treatment, Transplant_Treatment, Rep, Scaled_MassLoss)

# Average replicates within each year 
decomp_2021_avg <- decomp_2021_scaled %>%
  group_by(Site, CageID, Population, Trophic_Treatment, Transplant_Treatment) %>%
  summarise(Scaled_MassLoss_2021 = mean(Scaled_MassLoss, na.rm = TRUE), .groups = 'drop')

decomp_2023_avg <- decomp_2023_scaled %>%
  group_by(Site, CageID, Population, Trophic_Treatment, Transplant_Treatment) %>%
  summarise(Scaled_MassLoss_2023 = mean(Scaled_MassLoss, na.rm = TRUE), .groups = 'drop')


# Combine the 2021 and 2023 scaled data
combined_scaled_data <- decomp_2021_avg %>%
  left_join(decomp_2023_avg, by = c("Site", "CageID", "Population", "Trophic_Treatment", "Transplant_Treatment")) %>%
  mutate(Scaled_Difference = Scaled_MassLoss_2023 - Scaled_MassLoss_2021) %>% 
  filter(!is.na(Scaled_Difference) & !is.nan(Scaled_Difference))


# Plotting the change in decomposition rate across transplant treatments
ggplot(combined_scaled_data, aes(x = Population, y = Scaled_Difference, fill = Transplant_Treatment, color = Transplant_Treatment)) +
  geom_boxplot(alpha = 0.3, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.75), size = 1.1) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_classic() +
  coord_flip() +
  labs(x = "Population", y = "Difference in Mass Loss Rate (g/month)") +
  ggtitle("Scaled Difference in Decomposition Rate (Treatment - Baseline)")
  
# Plotting the change in decomposition rate across trophic treatments
ggplot(combined_scaled_data, aes(x = Population, y = Scaled_Difference, fill = Trophic_Treatment, color = Trophic_Treatment)) +
  geom_boxplot(alpha = 0.3, position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.75), size = 1.1) +
  scale_fill_brewer(palette = "Paired") +
  scale_color_brewer(palette = "Paired") +
  theme_classic() +
  coord_flip() +
  labs(x = "Population", y = "Difference in Mass Loss Rate (g/month)") +
  ggtitle("Scaled Difference in Decomposition Rate (Treatment - Baseline)")



### Regional variation: 2021 ----

# Visualize scaled mass loss
ggplot(decomp_2021_scaled, aes(x = Site, y = Scaled_MassLoss, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(width = 0.2, size = 1.1, alpha = 0.3) +
  theme_classic() +
  coord_flip() +
  labs(x = "Site", y = "Scaled Mass Loss (g/month)") +
  ggtitle("Scaled Decomposition Rates Across Sites")

# Calculate variability within and across sites
variability_within_sites <- decomp_2021_scaled %>%
  group_by(Site) %>%
  summarise(Within_Site_SD = sd(Scaled_MassLoss, na.rm = TRUE))

overall_sd <- sd(decomp_2021_scaled$Scaled_MassLoss, na.rm = TRUE)

print(variability_within_sites)
print(overall_sd)

# Field-visualization

# Calculate sample size for each site
sample_sizes <- decomp_2021 %>%
  group_by(Site) %>%
  summarise(Sample_Size = n())

variability_within_sites <- variability_within_sites %>%
  mutate(Variability_Type = "Within Site")

overall_variability <- data.frame(
  Site = "Overall",
  Within_Site_SD = overall_sd,
  Variability_Type = "Across Sites"
)

# Combine the data frames
variability_data <- bind_rows(variability_within_sites, overall_variability)

# Merge with sample sizes
variability_data <- variability_data %>%
  left_join(sample_sizes, by = "Site") %>%
  mutate(Sample_Size = ifelse(is.na(Sample_Size), sum(sample_sizes$Sample_Size), Sample_Size))

# Reorder the Site factor levels
variability_data <- variability_data %>%
  mutate(Site = fct_reorder(Site, Within_Site_SD, .desc = TRUE))

# Visualize variability
ggplot(variability_data, aes(x = Site, y = Within_Site_SD, fill = Variability_Type)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
  geom_text(aes(label = Sample_Size), 
            position = position_dodge(width = 0.9), 
            vjust = 0.1, hjust = -0.2) +
  theme_classic() +
  labs(x = "Site", y = "Standard Deviation of Scaled Mass Loss (g/month)") +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  ggtitle("Variability of Scaled Mass Loss Within and Across Sites")

