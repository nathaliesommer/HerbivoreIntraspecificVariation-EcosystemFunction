# SoilMicros

library(tidyverse)

SoilMicros <- read.csv("Data/Elements/SoilMicros.csv")

SoilMicros <- 
  SoilMicros %>%
  mutate(Latitude = case_when(
    (Site == "FN") ~ "South",
    (Site == "YF") ~ "South",
    (Site == "SC") ~ "South",
    (Site == "HF") ~ "South",
    (Site == "SP") ~ "North",
    (Site == "DC") ~ "North",
    (Site == "MC") ~ "North",
    (Site == "UP") ~ "North"))

## ALL Site-Level Total 
SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Total") %>% 
  filter(!Element == "B") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

## ALL Site-Level Available Micros
SoilMicros %>% 
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(!Element == "B") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  #geom_jitter(width = .2, color = "dark gray", size = 1.1) +
  #theme_classic() +
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available Soil Micronutrients") +
  theme_classic()


## Individual elements
SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Al") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()


SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Ca") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Fe") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "K") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Mg") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Na") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "P") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "S") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

SoilMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Digestion == "Available") %>% 
  filter(Element == "Zn") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Total Soil Micronutrients") +
  theme_classic()

