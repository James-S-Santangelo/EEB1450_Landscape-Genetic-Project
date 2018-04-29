# Load packages
library(tidyr)
library(dplyr)

# Load habitat type data
Structure_data <- "data-raw/Indiv_Habitat.csv"
Structure_data <- read.csv(Structure_data, header = TRUE)

# Clean and filter
names(Structure_data)[names(Structure_data) == "Population"] <- "City_pop"
names(Structure_data)[names(Structure_data) == "Enviro"] <- "Habitat"
Structure_data$City_pop <- as.character(Structure_data$City_pop)

# Create Final dataset
Structure_data <- Structure_data %>% 
  separate(City_pop, sep = "(?<=[A-Za-z])(?=[0-9])", into = c("city_code", "pop"), remove = TRUE) %>%# Split Population code into City ID and Pop columns
  mutate(PlantID = paste(city_code, pop, Individual, sep = ".")) %>% # Add Plant ID
  select(PlantID, everything()) %>% # Move Plant ID to first column
  filter(City == "Fergus" | City == "Acton") %>% # Only Acton and Fergus
  filter(pop != 1) %>%
  select(PlantID, City, pop, Habitat) %>% # Select Plant ID and habitat columns
  as.data.frame()

# Write cleanes up dataset to CSV
write.csv(Structure_data, "data-clean/Structure-data.csv", row.names = FALSE)


