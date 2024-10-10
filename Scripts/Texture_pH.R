## Code by NR Sommer

library(tidyverse)

### Initial Conditions (2021) ####

### Field-Level ----

texture_ph <- read.csv("Data/Texture_pH/2021_Texture_pH.csv")

texture_ph <- texture_ph %>% 
  rename(Silt = X.Silt, Sand = X.Sand, Clay = X.Clay)

texture_ph <- 
  texture_ph %>% 
  mutate(Latitude = case_when(
    (Field == "FN") ~ "South",
    (Field == "YF") ~ "South",
    (Field == "SC") ~ "South",
    (Field == "DC") ~ "North",
    (Field == "MC") ~ "North",
    (Field == "UP") ~ "North"))

# pH by field
texture_ph %>% 
  mutate(Field = factor(Field, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  ggplot(aes(y = pH, x = Field)) + 
  geom_boxplot(fill = "light gray") +  
  #geom_bar(position = "dodge", stat = "summary", fun = mean, fill = "gray") +
  geom_jitter(width = .2, color = "black", size = 1.1) +
  coord_flip() + 
  theme_classic() + 
  # geom_vline(xintercept = 3.5, linetype="dotted") + 
  labs(x="Site", title = "pH by Site")

# Texture by field
texture_ph %>%
  pivot_longer(cols = c('Sand', 'Silt', 'Clay'), 
               names_to = 'GrainType',
               values_to = 'Percent') %>%
  mutate(Field = factor(Field, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>% 
  mutate(GrainType = factor(GrainType, level = c("Clay", "Silt", "Sand"))) %>%
  ggplot(aes(y = Percent, x = Field, fill = GrainType)) + 
  geom_bar(position = "stack", stat = "summary", fun = mean, alpha = 0.9) +
  coord_flip() + 
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Site", title = "Texture by Site") +
  scale_fill_brewer(palette="Greys")
  #geom_vline(xintercept = 3.5, linetype="dotted")

View(texture_ph %>%
  pivot_longer(cols = c('Sand', 'Silt', 'Clay'), 
               names_to = 'GrainType',
               values_to = 'Percent') %>%
  group_by(Field, GrainType) %>%
  mutate(GrainTypeAvg = mean(Percent)))
