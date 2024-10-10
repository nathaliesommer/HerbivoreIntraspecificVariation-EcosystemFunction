# Plant Micros

library(tidyverse)

PlantMicros <- read.csv("Data/Elements/PlantMicros.csv")

PlantMicros <- 
  PlantMicros %>%
  mutate(Latitude = case_when(
    (Site == "FN") ~ "South",
    (Site == "YF") ~ "South",
    (Site == "SC") ~ "South",
    (Site == "HF") ~ "South",
    (Site == "SP") ~ "North",
    (Site == "DC") ~ "North",
    (Site == "MC") ~ "North",
    (Site == "UP") ~ "North"))

### SORU ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Aluminum ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Al") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Boron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "B") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()


#### Calcium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Ca") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Iron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Fe") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Potassium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "K") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Magnesium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Mg") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Manganese ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Mn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Sodium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Na") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Phosphorous ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "P") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Sulfer ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "S") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

#### Zinc ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "SORU") %>%
  filter(Element == "Zn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available SORU Micronutrients") +
  theme_classic()

### POPRC ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Aluminum ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Al") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Boron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "B") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()


#### Calcium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Ca") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Iron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Fe") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Potassium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "K") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Magnesium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Mg") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Manganese ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Mn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Sodium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Na") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Phosphorous ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "P") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Sulfer ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "S") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

#### Zinc ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "POPRC") %>%
  filter(Element == "Zn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available POPRC Micronutrients") +
  theme_classic()

### MISC ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Aluminum ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Al") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Boron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "B") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()


#### Calcium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Ca") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Iron ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Fe") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Potassium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "K") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Magnesium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Mg") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Manganese ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Mn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Sodium ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Na") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Phosphorous ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "P") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Sulfer ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "S") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()

#### Zinc ----
PlantMicros %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF","SP", "DC", "MC", "UP"))) %>%
  filter(Plant == "MISC") %>%
  filter(Element == "Zn") %>% 
  ggplot(aes(y=Concentration_mg.kg, x=Site)) +
  geom_boxplot(alpha=0.3, fill = "dark gray") +  
  coord_flip() +
  facet_wrap(~Element, nrow = 5) +
  labs(x="Site", y = "mg/kg", title = "Site-Level Available MISC Micronutrients") +
  theme_classic()



## Make a graph that compares the same element across plant types and across sites 
