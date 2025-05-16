# Texture, pH, bulk density

# Texture, pH, and bulk density were assessed at the site-level in Year 1, not in mesocosms

library(dplyr)

# Data
bulk_density <- read.csv("Data/BulkDensity/2021_BulkDensity.csv")
texture_ph <- read.csv("Data/Texture_pH/2021_Texture_pH.csv")

# Calculate bulk density
bulk_density_processed <- bulk_density %>% 
  mutate(SoilMass = Weight_SoilTray - Weight_Tray) %>% 
  mutate(SoilMassAdj = SoilMass - (Materials_volume_mL_end - Materials_volume_mL_start)) %>% 
  mutate(SoilBulkDensity = (SoilMassAdj / 308.89)) %>% # units = g/cm^3
  filter(Site %in% c("FN", "UP", "YF", "SC", "DC"))

# Process texture and pH data
texture_ph_processed <- texture_ph %>%
  filter(Field %in% c("FN", "UP", "YF", "SC", "DC"))

# Calculate summary statistics
bulk_density_summary <- bulk_density_processed %>%
  group_by(Site) %>%
  summarise(
    BD_mean = mean(SoilBulkDensity, na.rm = TRUE),
    BD_var = var(SoilBulkDensity, na.rm = TRUE)
  )

texture_ph_summary <- texture_ph_processed %>%
  group_by(Field) %>%
  summarise(
    pH_mean = mean(pH, na.rm = TRUE),
    pH_var = var(pH, na.rm = TRUE),
    Sand_mean = mean(X.Sand, na.rm = TRUE),
    Sand_var = var(X.Sand, na.rm = TRUE),
    Silt_mean = mean(X.Silt, na.rm = TRUE),
    Silt_var = var(X.Silt, na.rm = TRUE),
    Clay_mean = mean(X.Clay, na.rm = TRUE),
    Clay_var = var(X.Clay, na.rm = TRUE)
  )

print("Bulk Density Summary (g/cm³):")
print(bulk_density_summary)
print("\nTexture and pH Summary:")
print(texture_ph_summary)

# Calculate aggregate statistics across all sites
aggregate_stats <- list(
  "Bulk Density (g/cm³)" = list(
    mean = mean(bulk_density_processed$SoilBulkDensity, na.rm = TRUE),
    sd = sd(bulk_density_processed$SoilBulkDensity, na.rm = TRUE)
  ),
  "pH" = list(
    mean = mean(texture_ph_processed$pH, na.rm = TRUE),
    sd = sd(texture_ph_processed$pH, na.rm = TRUE)
  ),
  "Sand (%)" = list(
    mean = mean(texture_ph_processed$X.Sand, na.rm = TRUE),
    sd = sd(texture_ph_processed$X.Sand, na.rm = TRUE)
  ),
  "Silt (%)" = list(
    mean = mean(texture_ph_processed$X.Silt, na.rm = TRUE),
    sd = sd(texture_ph_processed$X.Silt, na.rm = TRUE)
  ),
  "Clay (%)" = list(
    mean = mean(texture_ph_processed$X.Clay, na.rm = TRUE),
    sd = sd(texture_ph_processed$X.Clay, na.rm = TRUE)
  )
)

# Aggregate statistics
for(var in names(aggregate_stats)) {
  cat(sprintf("%s: %.2f ± %.2f\n", 
              var, 
              aggregate_stats[[var]]$mean, 
              aggregate_stats[[var]]$sd))
}

# Plots
library(ggplot2)
library(tidyr)
library(paletteer)

bd_ph_plot_data <- bind_rows(
  bulk_density_processed %>%
    dplyr::select(Site, SoilBulkDensity) %>%
    rename(value = SoilBulkDensity) %>%
    mutate(variable = "Bulk Density (g/cm³)"),
  texture_ph_processed %>%
    dplyr::select(Field, pH) %>%
    rename(Site = Field, value = pH) %>%
    mutate(variable = "pH")
)

texture_plot_data <- texture_ph_processed %>%
  dplyr::select(Field, X.Sand, X.Silt, X.Clay) %>%
  rename(Site = Field) %>%
  pivot_longer(cols = c(X.Sand, X.Silt, X.Clay),
               names_to = "variable",
               values_to = "value") %>%
  mutate(variable = case_when(
    variable == "X.Sand" ~ "Sand",
    variable == "X.Silt" ~ "Silt",
    variable == "X.Clay" ~ "Clay"
  ))

# Set factor levels
texture_plot_data <- texture_plot_data %>%
  mutate(variable = factor(variable, levels = c("Sand", "Silt", "Clay")))

# pH and bulk density plot
ggplot(bd_ph_plot_data, aes(x = Site, y = value)) +
  geom_boxplot(aes(fill = Site)) +
  facet_wrap(~variable, scales = "free_y") +
  theme_bw() +
  labs(y = "Value", x = "Site") +
  theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 10)) +
  scale_fill_paletteer_d("NatParksPalettes::Yellowstone")


ggplot(texture_plot_data, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = Site)) +
  theme_bw() +
  labs(y = "Percent (%)", x = "Soil Fraction", fill = "Site") +
  scale_fill_paletteer_d("NatParksPalettes::Yellowstone")



