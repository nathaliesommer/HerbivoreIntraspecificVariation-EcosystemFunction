#### SIR calculation 
## Code by A. Polussa, A.M. Dobson, N.R. Sommer
# UPDATED Sep 17 2024

library(tidyverse)
library(dplyr)

# Function to calculate mean standard values
calculate_mean_standard <- function(stds) {
  stds %>% 
    rowwise() %>%
    mutate(meanStandard = mean(
      c_across(starts_with("std.value")), na.rm = TRUE
    )) %>%
    ungroup()
}

# Function to interpolate corrected standards
interpolate_corrected_standards <- function(stds, sir) {
  corr_std <- numeric()
  the.time <- numeric()
  the.slope <- numeric()
  v.num <- 1

  for (j in 2:nrow(stds)) {
    for (i in 1:(stds$irga.id[j] - stds$irga.id[j-1])) {
      time_diff_stds <- as.numeric(
        difftime(
          as.POSIXct(stds$std.start.time[j], format = "%m/%d/%Y %H:%M"),
          as.POSIXct(stds$std.start.time[j-1], format = "%m/%d/%Y %H:%M"), 
          units = "min"
        )
      )
      time_diff_sample <- as.numeric(
        difftime(
          as.POSIXct(sir$time.irga[stds$irga.id[j-1] + i - 1], 
                     format = "%m/%d/%Y %H:%M"),
          as.POSIXct(sir$time.irga[stds$irga.id[j-1]], 
                     format = "%m/%d/%Y %H:%M"), 
          units = "min"
        )
      )

      corr_std[v.num] <- as.numeric(stds$meanStandard[j-1]) +  
        ((as.numeric(stds$meanStandard[j]) - 
          as.numeric(stds$meanStandard[j-1])) / time_diff_stds) * 
        time_diff_sample

      the.time[v.num] <- as.numeric(
        difftime(
          as.POSIXct(sir$time.irga[stds$irga.id[j-1] + i - 1], format = "%m/%d/%Y %H:%M"),
          as.POSIXct(stds$std.end.time[j-1], format = "%m/%d/%Y %H:%M"),
          units = "min"
        )
      )

      the.slope[v.num] <- (
        (as.numeric(stds$meanStandard[j]) - as.numeric(stds$meanStandard[j-1])) / 
        time_diff_stds
      )
      v.num <- v.num + 1
    }
  }

  list(corr_std = corr_std, the.time = the.time, the.slope = the.slope)
}

# Function to calculate SIR
calculate_sir <- function(sir, gwc) {
  stds <- sir %>% filter(!is.na(std.value.1)) %>% calculate_mean_standard()
  interpolated <- interpolate_corrected_standards(stds, sir)

  sir <- sir[1:(nrow(sir) - 1),] # remove last row

  sir <- sir %>%
    mutate(
      correctedStandard = as.numeric(interpolated$corr_std),
      the.time = interpolated$the.time,
      the.slope = interpolated$the.slope
    ) %>%
    dplyr::select(-starts_with("std.value"), -std.start.time, -std.end.time)

  sir_calc <- sir %>%
    mutate(
      incubationTime = as.numeric(
        difftime(
          as.POSIXct(time.irga, format = "%m/%d/%Y %H:%M"),
          as.POSIXct(time.flush, format = "%m/%d/%Y %H:%M"), 
          units = "hours"
        )
      ),
      dilutionFactor = ((5 * times.sampled) / (57.15 - soil.volume)) + 1,
      measuredCO2 = irga.integral * (standard.co2 / correctedStandard),
      concentrationCO2 = measuredCO2 * dilutionFactor,
      volumeCO2 = concentrationCO2 * ((57.15 - soil.volume) / 1000),
      molesCO2 = (volumeCO2 / 22.414) * 273.15 / 293.15,
      CO2C = molesCO2 * 12.011,
      CO2CperHour = CO2C / incubationTime
    ) %>%
    dplyr::select(
      irga.id, unique.id, replicate, soil.volume, actual.fresh.mass, 
      standard.co2, correctedStandard, the.time, the.slope, irga.integral, 
      incubationTime, CO2C, CO2CperHour
    )

  sir_calc_normalized <- sir_calc %>%
    left_join(gwc, by = "unique.id") %>%
    mutate(
      DryMass = actual.fresh.mass * (1 - moistureFraction),
      CO2CperHourperg = CO2CperHour / DryMass
    ) %>%
    dplyr::select(unique.id, replicate, CO2CperHourperg)

  sir_calc_normalized %>%
    group_by(unique.id) %>%
    summarize(CO2CperHourperg = mean(CO2CperHourperg, na.rm = TRUE))
}

# Main script
sir_files_2021 <- c(
  "Data/SIR/2021_SIR/soil_sir_11.21.csv", 
  "Data/SIR/2021_SIR/soil_sir_11.27.csv", 
  "Data/SIR/2021_SIR/soil_sir_11.28.csv", 
  "Data/SIR/2021_SIR/soil_sir_11.29.csv", 
  "Data/SIR/2021_SIR/soil_sir_12.27.csv", 
  "Data/SIR/2021_SIR/soil_sir_12.29.csv", 
  "Data/SIR/2021_SIR/soil_sir_1.1.csv"
)
sir_files_2023 <- c(
  "Data/SIR/2023_SIR/soil_sir_11.08.csv", 
  "Data/SIR/2023_SIR/soil_sir_11.14.csv", 
  "Data/SIR/2023_SIR/soil_sir_11.17.csv", 
  "Data/SIR/2023_SIR/soil_sir_11.19.csv"
)
gwc_file_2021 <- "Data/SIR/2021_SIR/soilGWC_2021.csv"
gwc_file_2023 <- "Data/SIR/2023_SIR/soilGWC_2023.csv"

process_sir_data <- function(sir_files, gwc_file) {
  all_sir_data <- list()
  
  for (sir_file in sir_files) {
    sir <- read_csv(sir_file)
    gwc <- read_csv(gwc_file)
    sir_data <- calculate_sir(sir, gwc)
    all_sir_data <- bind_rows(all_sir_data, sir_data)
  }
  
  all_sir_data
}

# Process 2021 data
sir_data_2021 <- process_sir_data(sir_files_2021, gwc_file_2021)

# Process 2023 data
sir_data_2023 <- process_sir_data(sir_files_2023, gwc_file_2023)

# Combine with additional data for 2021
prep_data_2021 <- read.csv("Data/SIR/2021_SIR/SIR Prep_2021.csv") %>%
  dplyr::select(Sample_ID, unique.id)

c_2021 <- left_join(sir_data_2021, prep_data_2021, by = "unique.id") %>%
  mutate(Year = 2021) %>% 
  drop_na()
  

# Combine with additional data for 2023
prep_data_2023 <- read.csv("Data/SIR/2023_SIR/SIR Prep_2023.csv") %>%
  dplyr::select(Sample_ID, unique.id)

c_2023 <- left_join(sir_data_2023, prep_data_2023, by = "unique.id") %>%
  mutate(Year = 2023) %>% 
  drop_na()

# Calculate the average CO2CperHourperg for each Sample_ID within each year
SIR_avg_2021 <- c_2021 %>%
  group_by(Sample_ID, Year) %>%
  summarize(CO2CperHourperg = mean(CO2CperHourperg, na.rm = TRUE))

SIR_avg_2023 <- c_2023 %>%
  group_by(Sample_ID, Year) %>%
  summarize(CO2CperHourperg = mean(CO2CperHourperg, na.rm = TRUE))

# Combine the 2021 and 2023 data frames by Sample_ID and Year
SIR_final_data <- bind_rows(SIR_avg_2021, SIR_avg_2023) %>%
  filter(!is.na(Sample_ID))

# Apply additional transformations
SIR_final_data <- SIR_final_data %>%
  mutate(Site = case_when(
    startsWith(as.character(Sample_ID), "FNH") ~ "FN",
    startsWith(as.character(Sample_ID), "FNS") ~ "YF",
    startsWith(as.character(Sample_ID), "FNN") ~ "UP",
    startsWith(as.character(Sample_ID), "MCH") ~ "MC",
    startsWith(as.character(Sample_ID), "MCN") ~ "UP",
    startsWith(as.character(Sample_ID), "MCS") ~ "YF",
    startsWith(as.character(Sample_ID), "YFH") ~ "YF",
    startsWith(as.character(Sample_ID), "YFN") ~ "UP",
    startsWith(as.character(Sample_ID), "SCN") ~ "UP",
    startsWith(as.character(Sample_ID), "SCH") ~ "SC",
    startsWith(as.character(Sample_ID), "SCS") ~ "YF",
    startsWith(as.character(Sample_ID), "UPH") ~ "UP",
    startsWith(as.character(Sample_ID), "UPS") ~ "YF",
    startsWith(as.character(Sample_ID), "DCH") ~ "DC",
    startsWith(as.character(Sample_ID), "DCN") ~ "UP",
    startsWith(as.character(Sample_ID), "DCS") ~ "YF"
  )) %>%
  mutate(Trophic_Treatment = case_when(
    substr(Sample_ID, 5, 5) == "V" ~ "Vegetation",
    substr(Sample_ID, 5, 5) == "H" ~ "Herbivore",
    substr(Sample_ID, 5, 5) == "P" ~ "Predator"
  )) %>%
  mutate(Population = substr(Sample_ID, 1, 2))


### Permute UPH_H7 and YFN_H2 2023

# Create new rows with specified information
new_sir_rows <- tibble(
  Sample_ID = c("UPH_H7", "YFN_H2"),
  Year = c(2023, 2023),
  Population = c("UP", "YF"),
  Site = c("UP", "UP"),
  Trophic_Treatment = c("Herbivore", "Herbivore")
)

# Function to calculate average CO2CperHourperg based on matching conditions
average_sir_rates <- function(data, new_row) {
  # Filter existing data to find matching rows
  matching_rows <- data %>%
    filter(
      Year == new_row$Year,
      Population == new_row$Population,
      Site == new_row$Site,
      Trophic_Treatment == new_row$Trophic_Treatment
    )
  
  # Calculate the average CO2CperHourperg if matching rows are found
  if (nrow(matching_rows) > 0) {
    CO2CperHourperg <- mean(matching_rows$CO2CperHourperg, na.rm = TRUE)
  } else {
    # If no matching rows are found, set CO2CperHourperg to NA
    CO2CperHourperg <- NA
  }
  
  return(CO2CperHourperg)
}

# Apply the average calculation function to each new row
new_sir_rows <- new_sir_rows %>%
  rowwise() %>%
  mutate(
    CO2CperHourperg = average_sir_rates(SIR_final_data, list(
      Year = Year,
      Population = Population,
      Site = Site,
      Trophic_Treatment = Trophic_Treatment
    ))
  ) %>%
  ungroup()

# Combine the new rows with the existing data
SIR_final_data <- bind_rows(SIR_final_data, new_sir_rows)




