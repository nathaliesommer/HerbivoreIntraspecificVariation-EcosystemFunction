# Load necessary libraries
library(dplyr)
library(lubridate)

# Define the directory containing the data files
data_dir <- "Data/Temperatures"

# List all files in the directory
all_files <- list.files(data_dir, full.names = TRUE)

# Function to read and standardize data with error handling
read_and_standardize <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path)
    
    # Standardize column names
    colnames(data) <- tolower(colnames(data))
    colnames(data) <- gsub("\\.", "_", colnames(data))
    
    # Handle different temperature column names
    temp_col <- grep("temperaturec|temperature_c|temp_c|temp", colnames(data), value = TRUE)
    if (length(temp_col) == 0) {
      stop("No temperature column found.")
    }
    
    # Rename temperature column to a standard name
    colnames(data)[colnames(data) == temp_col] <- "temperaturec"
    
    # Define the required columns
    required_columns <- c("event", "date_time", "temperaturec", "field", "pop", "cage", "treatment")
    
    # Drop any extra columns
    data <- data[ , required_columns, drop = FALSE]
    
    # Convert date_time to Date format using lubridate for flexibility
    data$date <- as.Date(parse_date_time(data$date_time, orders = c("Ymd HMS", "mdy HM", "mdy HMS", "ymd HM", "ymd HMS", "mdy HMS")))
    
    return(list(data = data, temp_col = "temperaturec"))
  }, error = function(e) {
    message("Error reading or processing file: ", file_path, " - ", e$message)
    return(NULL)
  })
}

# Function to combine data across all loggers for each field and year
combine_data_across_loggers <- function(logger_groups, date_ranges) {
  combined_data <- list()
  
  for (logger_id in names(logger_groups)) {
    logger_files <- logger_groups[[logger_id]]
    
    # Combine data from all files for this logger
    logger_data <- do.call(rbind, lapply(logger_files, function(file_path) {
      result <- read_and_standardize(file_path)
      if (is.null(result)) return(NULL)
      return(result$data)
    }))
    
    # Add logger data to combined data
    combined_data[[logger_id]] <- logger_data
  }
  
  # Bind all logger data into a single data frame
  combined_data <- bind_rows(combined_data)
  
  # Print unique fields and years for debugging
  print(unique(combined_data$field))
  print(unique(year(combined_data$date)))
  
  # Filter out data for the year 2024
  combined_data <- combined_data %>% filter(year(date) != 2024)
  
  # Filter data by date ranges and add year column
  aggregated_data <- list()
  for (year in names(date_ranges)) {
    start_date <- as.Date(date_ranges[[year]][1])
    end_date <- as.Date(date_ranges[[year]][2])
    
    yearly_data <- combined_data %>%
      filter(date >= start_date & date <= end_date) %>%
      mutate(Year = year)
    
    if (nrow(yearly_data) > 0) {
      aggregated_data[[year]] <- yearly_data
    }
  }
  
  return(bind_rows(aggregated_data))
}

# Function to calculate metrics
calculate_metrics <- function(aggregated_data) {
  # Calculate daily max for each field and date
  daily_max <- aggregated_data %>%
    group_by(field, Year, date) %>%
    summarise(
      Daily_Max = if (all(is.na(temperaturec))) NA else max(temperaturec, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Calculate metrics for each field and year
  metrics <- daily_max %>%
    group_by(field, Year) %>%
    summarise(
      MeanDailyMax = mean(Daily_Max, na.rm = TRUE),
      AverageTemp = mean(aggregated_data$temperaturec[aggregated_data$field == first(field) & aggregated_data$Year == first(Year)], na.rm = TRUE),
      Coefficient_of_Variation = sd(aggregated_data$temperaturec[aggregated_data$field == first(field) & aggregated_data$Year == first(Year)], na.rm = TRUE) / 
        mean(aggregated_data$temperaturec[aggregated_data$field == first(field) & aggregated_data$Year == first(Year)], na.rm = TRUE)
    )
  
  return(metrics)
}

# Extract logger identifier from file names
extract_logger_id <- function(file_name) {
  sub("^(S_[^_]+_[^_]+_[^_]+)_.*\\.csv$", "\\1", basename(file_name))
}

# Define date ranges for each year
date_ranges <- list(
  "2021" = c("2021-06-15", "2021-10-01"),
  "2022" = c("2022-06-15", "2022-10-01"),
  "2023" = c("2023-06-15", "2023-10-01")
)

#### 60cm LOGGERS -------

# Filter files based on naming structure to include all files with '60'
selected_files <- grep("S_.*_60_.*\\.csv$", all_files, value = TRUE)

# Group files by logger identifier
logger_groups <- split(selected_files, sapply(selected_files, extract_logger_id))

# Main processing
combined_data <- combine_data_across_loggers(logger_groups, date_ranges)
final_output_60 <- calculate_metrics(combined_data)

# Print the final output
print("Final Output:")
print(final_output_60)

# write.csv(final_output_60, "TemperatureMetaData60.csv")

### ALL LOGGERS -------

# Main processing for all instances (not filtering by height)
all_files <- list.files(data_dir, full.names = TRUE)  # Re-list all files without filtering
logger_groups_all <- split(all_files, sapply(all_files, extract_logger_id))

combined_data_all <- combine_data_across_loggers(logger_groups_all, date_ranges)

final_output_all <- calculate_metrics(combined_data_all)

# Print the final output for all instances
print("Final Output for all instances:")
print(final_output_all)

#### HOME-SITE CAGES and NO PREDATORS (all heights) ----

selected_files <- grep("S_..H_[^P].*\\.csv$", all_files, value = TRUE)
exceptions <- file.path(data_dir, c("S_FNH_P1_20_Fall2021.csv", "S_FNH_P1_60_Fall2021.csv")) # only two loggers correctly calibrated in Fall 2021 at FN
selected_files <- c(selected_files, exceptions)


# Group files by logger identifier
logger_groups <- split(selected_files, sapply(selected_files, extract_logger_id))

# Main processing
combined_data <- combine_data_across_loggers(logger_groups, date_ranges)
final_output_cleaned <- calculate_metrics(combined_data)

# Print the final output
print("Final Output:")
print(final_output_cleaned)
# write.csv(final_output_cleaned, "Data/Temperatures/TemperatureMetaData.csv")

#### AVERAGED ACROSS YEARS 2022 AND 2023 -------
# **Update with FNH temperature data**
# Home-site cages, no predators, all heights, temperature metrics across years 2022 and 2023

# Filter data for the years 2022 and 2023
combined_data_2022_2023 <- combined_data %>%
  filter(Year %in% c("2022", "2023"))

# Calculate daily max for each logger, field, and date
daily_max_2022_2023 <- combined_data_2022_2023 %>%
  group_by(field, date) %>%
  summarise(
    Daily_Max = if (all(is.na(temperaturec))) NA else max(temperaturec, na.rm = TRUE),
    .groups = 'drop'
  )

# Calculate metrics for the combined years
temp_metrics <- combined_data_2022_2023 %>%
  group_by(field) %>%
  summarise(
    AverageTemp = mean(temperaturec, na.rm = TRUE),  # True average temperature
    MeanDailyMax = mean(daily_max_2022_2023$Daily_Max[daily_max_2022_2023$field == first(field)], na.rm = TRUE),  # Average of daily max temperatures
    Coefficient_of_Variation = sd(temperaturec, na.rm = TRUE) / mean(temperaturec, na.rm = TRUE)  # Coefficient of variation
  )

# Print the metrics for the combined years
print("Metrics for combined years 2022 and 2023:")
print(temp_metrics)
temp_metrics <- temp_metrics %>% 
  rename(Site = field) 
