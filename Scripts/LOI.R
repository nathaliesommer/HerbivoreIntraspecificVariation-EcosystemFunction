## Loss-on-ignition

# Load necessary library
library(dplyr)

# Read in the CSV file
data <- read.csv("Data/LOI/2021_LOI.csv")

# Process the data
SOM_data <- data %>%
  # Drop columns where Comments = "sample error"
  filter(Comments != "sample error") %>%
  
  # Drop rows where Site is NA
  filter(!is.na(Site)) %>%
  
  # Take the average of "Adjusted_SOM_Percent" for each combination of Sample_ID
  group_by(Sample_ID) %>%
  summarise(average_SOM = mean(Adjusted_SOM_Percent, na.rm = TRUE)) %>%
  
  # Drop the first two characters of the Sample_ID
  mutate(Sample_ID = substr(Sample_ID, 3, nchar(Sample_ID))) %>%
  
  # Drop cases where the first and second characters = "HF" or "MC"
  filter(!grepl("^(HF|MC)", Sample_ID)) %>%
  
  # Drop cases where the fifth character = "P"
  filter(substr(Sample_ID, 5, 5) != "P")

# View the processed data
print(SOM_data)

