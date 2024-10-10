### Bulk Density -
library(tidyverse)

BulkDensity <- read.csv("Data/BulkDensity/2021_BulkDensity.csv")

BulkDensity <- BulkDensity %>% 
  mutate(SoilMass = Weight_SoilTray-Weight_Tray) %>% 
  mutate(SoilMassAdj = SoilMass - (Materials_volume_mL_end - Materials_volume_mL_start)) %>%
  mutate(SoilBulkDensity = 
           (SoilMassAdj / 308.89)) # units = g/cm^3, where 308.89 is the volume of the bulk density corer

# PLOTS

BulkDensity <- 
  BulkDensity %>% 
  mutate(Latitude = case_when(
    (Site == "FN") ~ "South",
    (Site == "YF") ~ "South",
    (Site == "SC") ~ "South",
    (Site == "DC") ~ "North",
    (Site == "MC") ~ "North",
    (Site == "UP") ~ "North"))

BulkDensity %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  ggplot(aes(y=SoilBulkDensity, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  geom_jitter(width = .2, color = "black", size = 1.1) +
  theme_classic() +
  coord_flip() +
  #geom_vline(xintercept = 3.5, linetype = "dotted") + 
  labs(x="Site", y = "Bulk Density (g/cm3)")
