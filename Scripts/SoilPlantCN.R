# CarbonNitrogen ----

library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(paletteer)

CNdata <- read.csv("Data/Elements/CN_Tidy.csv")

CNdata <- CNdata %>% 
  mutate(
    PercentN = ifelse(Note %in% c("None", "Defoliated"), 0, PercentN),
    PercentC = ifelse(Note %in% c("None", "Defoliated"), 0, PercentC)
  )

CNdata <- CNdata %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))